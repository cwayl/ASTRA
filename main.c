#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <string.h>


#define ARRAYSIZE(a) (sizeof(a) / sizeof((a)[0]))
#define NUM_ROWS(a) ARRAYSIZE(a)
#define NUM_COLS(a) ARRAYSIZE((a)[0])
#define MAX_ITERATIONS 100
#define TOLERANCE 1E-10
#define NUM_ROWS_Q 3 // Number of rows in Q
#define NUM_COLS_Q 3 // Number of columns in Q
#define NUM_ROWS_RHO_G 3 // Number of rows in rho_G (position vector in perifocal frame)
#define NUM_COLS_RHO_G 1 // Number of columns in rho_G, which is 1 since it's a vector

// Function prototypes
double solveKeplersEquation(double meanAnomaly, double eccentricity);
double calculateTrueAnomaly(double eccentricAnomaly, double eccentricity);

int main() {

    int socket_desc;
    struct sockaddr_in server_addr;
    char server_message[2000];

    // Create socket:
    socket_desc = socket(AF_INET, SOCK_STREAM, 0);

    if(socket_desc < 0){
        printf("Unable to create socket\n");
        return -1;
    }
    printf("Socket created successfully\n");

    // Set port and IP the same as server-side:
    server_addr.sin_family = AF_INET;
    server_addr.sin_port = htons(4533);
    server_addr.sin_addr.s_addr = inet_addr("127.0.0.1");

    // Send connection request to server:
    if(connect(socket_desc, (struct sockaddr*)&server_addr, sizeof(server_addr)) < 0){
        printf("Unable to connect\n");
        return -1;
    }
    printf("Connected with server successfully\n");

    // Start Cooper Code

    // Define variables for user input
    double latitude, longitude, altitude;
    int epochYear;
    double epoch, inclination, raan, eccentricity, argumentOfPerigee, meanAnomaly, meanMotion;


    // Step 1: Prompt user for TLE input
    /*
    printf("Enter Latitude of Ground Station (degrees): ");
    scanf("%lf", &latitude);

    printf("Enter Longitude of Ground Station (degrees): ");
    scanf("%lf", &longitude);

    printf("Enter Altitude of Ground Station (km): ");
    scanf("%lf", &altitude);

    printf("Enter Epoch Year from TLE (YYYY): ");
    scanf("%d", &epochYear);

    printf("Enter Epoch from TLE (day of the year): ");
    scanf("%lf", &epoch);

    printf("Enter Inclination from TLE (degrees): ");
    scanf("%lf", &inclination);

    printf("Enter RAAN (Right Ascension of the Ascending Node) from TLE (degrees): ");
    scanf("%lf", &raan);

    printf("Enter Eccentricity from TLE: ");
    scanf("%lf", &eccentricity);

    printf("Enter Argument of Perigee from TLE (degrees): ");
    scanf("%lf", &argumentOfPerigee);

    printf("Enter Mean Anomaly from TLE (degrees): ");
    scanf("%lf", &meanAnomaly);

    printf("Enter Mean Motion from TLE (revolutions per day): ");
    scanf("%lf", &meanMotion);

     */


    latitude = 39.531254;
    longitude = -82.409393;
    altitude = 1.382;
    epochYear = 2024;
    epoch = 59.80121106;
    inclination = 97.4396;
    raan = 127.9486;
    eccentricity = 0.0012457;
    argumentOfPerigee = 306.0086;
    meanAnomaly = 53.9993;
    meanMotion = 15.19316519;


/*
    latitude = 39.531254;
    longitude = -82.409393;
    altitude = 1.382;
    epochYear = 2024;
    epoch = 62.43560187;
    inclination = 97.4392;
    raan = 130.5442;
    eccentricity = 0.0012317;
    argumentOfPerigee = 297.0116;
    meanAnomaly = 63.9850;
    meanMotion = 15.19378124;
*/


    // Constants
    const double mu = 398600; // Earth's gravitational parameter in km^3/s^2
    double a, h; // Semi-major axis, radius of orbit at perigee, and angular momentum
    double n; // Mean Motion in radians per second

    // Convert Mean Motion from revolutions per day to radians per second
    n = meanMotion * (2 * M_PI) / (24 * 3600); // Corrected variable

    // Step 2: Calculating orbital angular momentum
    a = pow(mu / pow(n, 2), 1.0 / 3.0);
    h = sqrt(mu * a * (1 - pow(eccentricity,2)));

    // start loop here

    while(1) {

        // Find UTC time
        // Structure ptr gets UTC time
        struct tm *ptr;
        time_t t;
        t = time(NULL);
        ptr = gmtime(&t);

        // Assign UTC values to individual variables
        double year;
        double month;
        double day;
        double UT;
        year = ptr->tm_year + 1900;
        month = ptr->tm_mon + 1;
        day = ptr->tm_mday;
        UT = (float) ptr->tm_hour + (float) ptr->tm_min / 60 + (float) ptr->tm_sec / 3600;

        /*
        // for testing
        month = 2;
        day = 29;
        UT = 20;
*/

        double dayOfYear = ptr->tm_yday + 1 + UT/24;
        // for testing
        dayOfYear = 60 + (double) 20/24;
        double delta_t = dayOfYear - epoch;
        double delta_MeanAnamoly = delta_t * meanMotion * 360;


        double newMeanAnomaly;
        newMeanAnomaly = meanAnomaly + delta_MeanAnamoly;

        newMeanAnomaly = fmod(newMeanAnomaly, 360);

        // Convert Mean Anomaly from degrees to radians
        double meanAnomalyRadians = newMeanAnomaly * M_PI / 180.0;

        // Solve Kepler's Equation for Eccentric Anomaly
        double eccentricAnomaly = solveKeplersEquation(meanAnomalyRadians, eccentricity);

        // Step 3: Calculate True Anomaly
        double trueAnomaly = calculateTrueAnomaly(eccentricAnomaly, eccentricity);

        // Step 4: Calculate the position vector in the perifocal frame
        double r_perifocal[NUM_ROWS_RHO_G][NUM_COLS_RHO_G];
        double theta = trueAnomaly; // True anomaly is the angle (theta)

        r_perifocal[0][0] = (h * h / mu) * (1 / (1 + eccentricity * cos(theta))) * cos(theta); // x-component
        r_perifocal[1][0] = (h * h / mu) * (1 / (1 + eccentricity * cos(theta))) * sin(theta); // y-component
        r_perifocal[2][0] = 0; // z-component

        // printf("Position vector in the perifocal frame: [%lf, %lf, %lf]\n", r_perifocal[0][0], r_perifocal[1][0], r_perifocal[2][0]);

        // Define the transformation matrix Q based on inclination, RAAN, and argument of perigee
        double raanRad = raan * M_PI / 180.0;
        double inclinationRad = inclination * M_PI / 180.0;
        double argumentOfPerigeeRad = argumentOfPerigee * M_PI / 180.0;

        /*
        double Q[NUM_ROWS_Q][NUM_COLS_Q] = {
                {cos(raanRad) * cos(argumentOfPerigeeRad) - sin(raanRad) * sin(argumentOfPerigeeRad) * cos(inclinationRad), -cos(raanRad) * sin(argumentOfPerigeeRad) - sin(raanRad) * cos(argumentOfPerigeeRad) * cos(inclinationRad), sin(raanRad) * sin(inclinationRad)},
                {sin(raanRad) * cos(argumentOfPerigeeRad) + cos(raanRad) * sin(argumentOfPerigeeRad) * cos(inclinationRad), -sin(raanRad) * sin(argumentOfPerigeeRad) + cos(raanRad) * cos(argumentOfPerigeeRad) * cos(inclinationRad), -cos(raanRad) * sin(inclinationRad)},
                {sin(argumentOfPerigeeRad) * sin(inclinationRad), cos(argumentOfPerigeeRad) * sin(inclinationRad), cos(inclinationRad)}
        };
        */
        // Step 5: Switching from the perifocal frame to the geocentric equatorial frame
        double R_1[NUM_ROWS_Q][NUM_COLS_Q] = {
                {cos(argumentOfPerigeeRad), sin(argumentOfPerigeeRad), 0},
                {-sin(argumentOfPerigeeRad), cos(argumentOfPerigeeRad), 0},
                {0, 0, 1},
        };
        double R_2[NUM_ROWS_Q][NUM_COLS_Q] = {
                {1, 0 ,0},
                {0, cos(inclinationRad), sin(inclinationRad)},
                {0, -sin(inclinationRad), cos(inclinationRad)}

        };
        double R_3[NUM_ROWS_Q][NUM_COLS_Q] = {
                {cos(raanRad), sin(raanRad), 0},
                {-sin(raanRad), cos(raanRad), 0},
                {0, 0, 1},
        };

        double R_1_2resultMatrix[NUM_ROWS_Q][NUM_COLS_Q]={
                {0, 0, 0},
                {0, 0, 0},
                {0, 0, 0}
        };

        double Q_Matrix[NUM_ROWS_Q][NUM_COLS_Q]={
                {0, 0, 0},
                {0, 0, 0},
                {0, 0, 0}
        };

        // Perform matrix multiplication R_1 * R_2
            for (int i = 0; i < NUM_ROWS_Q; i++) {
                for (int j = 0; j < NUM_COLS_Q; j++) {
                    R_1_2resultMatrix[i][j] = 0; // Initialize the element to 0
                    for (int k = 0; k < NUM_COLS_Q; k++) {
                        R_1_2resultMatrix[i][j] += R_1[i][k] * R_2[k][j];
                    }
                }
            }

        // Perform matrix multiplication (R_1 * R_2) * R_3
        for (int i = 0; i < NUM_ROWS_Q; i++) {
            for (int j = 0; j < NUM_COLS_Q; j++) {
                Q_Matrix[i][j] = 0; // Initialize the element to 0
                for (int k = 0; k < NUM_COLS_Q; k++) {
                    Q_Matrix[i][j] += R_1_2resultMatrix[i][k] * R_3[k][j];
                }
            }
        }

        double Q_Matrix_New[NUM_ROWS_Q][NUM_COLS_Q]= {
                {0, 0, 0},
                {0, 0, 0},
                {0, 0, 0}
        };

        // Transpose the Q_Matrix
        for (int i = 0; i < NUM_ROWS_Q; i++) {
            for (int j = 0; j < NUM_COLS_Q; j++) {
                Q_Matrix_New[j][i] = Q_Matrix[i][j];
            }
        }

        // Initialize the geocentric equatorial frame position vector
        double r_geocentric[NUM_ROWS_Q][NUM_COLS_RHO_G] = {
                {Q_Matrix_New[0][0] * r_perifocal[0][0] + Q_Matrix_New[0][1] * r_perifocal[1][0] +Q_Matrix_New[0][2] * r_perifocal[2][0]},
                {Q_Matrix_New[1][0] * r_perifocal[0][0] + Q_Matrix_New[1][1] * r_perifocal[1][0] +Q_Matrix_New[1][2] * r_perifocal[2][0]},
                {Q_Matrix_New[2][0] * r_perifocal[0][0] + Q_Matrix_New[2][1] * r_perifocal[1][0] +Q_Matrix_New[2][2] * r_perifocal[2][0]}
        };

        // printf("Position vector in the geocentric equatorial frame: [%lf, %lf, %lf]\n", r_geocentric[0][0], r_geocentric[1][0], r_geocentric[2][0]);

        // End cooper Code

        // Start CJ Code

        /* Step 6 - Calculate sidereal time
         * Find UTC time
         * Convert to Julian time
         * Calculate Greenwich sidereal time
         * Calculate local sidereal time */

        // Julian time.
        double J_O = 367 * year - floor((7 * floor((month + 9) / 12)) / 4) + floor((275 * month) / 9) + day + 1721013.5;
        double T_O = (J_O - 2451545) / 36525;

        // Greenwich sidereal time
        double greenwich_sidereal_time_0 =
                100.4606184 + 36000.77004 * T_O + 0.000387933 * pow(T_O, 2) - 2.583 * pow(10, -8) * pow(T_O, 3);
        double greenwich_sidereal_time = greenwich_sidereal_time_0 + 360.98564724 * UT / 24;

        // local sidereal time
        double sidereal_time = greenwich_sidereal_time + longitude;
        sidereal_time = fmod(sidereal_time, 360);

        /* Step 7 - Calculate Observer Position in Geocentric Equatorial Frame
         * Values of shape of earth are given and constant
         * Convert latitude and sidereal time to radians to be used in cos and sin
         * Compute observer position
         * */

        // flattening of the earth
        double f = 0.003353;
        // Equatorial radius of the earth in km
        double R_e = 6378;

        // Convert lat and sidereal to radians
        double latitudeRad = latitude / 180 * M_PI;
        sidereal_time = sidereal_time / 180 * M_PI;

        // Observer position
        double R[3][1];
        R[0][0] = (R_e / (sqrt(1 - (2 * f - pow(f, 2)) * pow(sin(latitudeRad), 2))) + altitude) * cos(latitudeRad) *
                  cos(sidereal_time);
        R[1][0] = (R_e / (sqrt(1 - (2 * f - pow(f, 2)) * pow(sin(latitudeRad), 2))) + altitude) * cos(latitudeRad) *
                  sin(sidereal_time);
        R[2][0] = ((R_e * pow((1 - f), 2)) / (sqrt(1 - (2 * f - pow(f, 2)) * pow(sin(latitudeRad), 2))) + altitude) *
                  sin(latitudeRad);

        /* Step 8 - Calculate relative position vector of sat to station
         * Subtract the Observer position from the satellite position
         * */

        double rho_G[3][1];
        rho_G[0][0] = r_geocentric[0][0] - R[0][0];
        rho_G[1][0] = r_geocentric[1][0] - R[1][0];
        rho_G[2][0] = r_geocentric[2][0] - R[2][0];

        /* Step 9 - Calculate azimuth and Elevation
         * Make rotation matrix
         * multiply rho by rotation matrix
         * find unit vector of relative position vector
         * convert unit vector to Azimuth + Elevation
         * */

        // Make the rotation matrix
        double Q_2[3][3];

        Q_2[0][0] = -sin(sidereal_time);
        Q_2[0][1] = cos(sidereal_time);
        Q_2[0][2] = 0;
        Q_2[1][0] = -sin(latitudeRad) * cos(sidereal_time);
        Q_2[1][1] = -sin(latitudeRad) * sin(sidereal_time);
        Q_2[1][2] = cos(latitudeRad);
        Q_2[2][0] = cos(latitudeRad) * cos(sidereal_time);
        Q_2[2][1] = cos(latitudeRad) * sin(sidereal_time);
        Q_2[2][2] = sin(latitudeRad);

        // Multiplies the arrays Q and rho_G to get new array rho_R
        double rho_R[NUM_ROWS(Q_2)][NUM_COLS(rho_G)];
        int i, j, k;

        // Initializing elements of matrix to 0.
        for (i = 0; i < NUM_ROWS(Q_2); ++i) {
            for (j = 0; j < NUM_COLS(rho_G); ++j) {
                rho_R[i][j] = 0;
            }
        }

        // Multiplying matrix array1 and array2 and storing in array.
        for (i = 0; i < NUM_ROWS(Q_2); ++i) {
            for (j = 0; j < NUM_COLS(rho_G); ++j) {
                for (k = 0; k < NUM_COLS(Q_2); ++k) {
                    rho_R[i][j] += Q_2[i][k] * rho_G[k][j];
                }
            }
        }

        // normalize rho_R
        double rho_size = sqrt(pow(rho_R[0][0], 2) + pow(rho_R[1][0], 2) + pow(rho_R[2][0], 2));
        rho_R[0][0] = rho_R[0][0] / rho_size;
        rho_R[1][0] = rho_R[1][0] / rho_size;
        rho_R[2][0] = rho_R[2][0] / rho_size;

        // convert unit vector to Azimuth + Elevation
        double Elevation = asin(rho_R[2][0]);
        double check = rho_R[0][0] / cos(Elevation);
        double Azimuth = acos(rho_R[1][0] / cos(Elevation));
        if (check < 0) {
            Azimuth = 2 * M_PI - Azimuth;
        }

        // Convert from radians to degrees
        double elevationRad = Elevation * 180 / M_PI;
        double azimuthRad = Azimuth * 180 / M_PI;

        /* Step 10 - Check
         * Check z height
         * */

        if(rho_R[0][0] < 0) {
            // printf("Move to Home\n");
        } else {
            // printf("Move to Coordinates\n");
        }

        printf("Azimuth: %f\n", azimuthRad);
        printf("Elevation: %f\n", elevationRad);

        int azimuthInt = floor(azimuthRad);
        int elevationInt = floor(elevationRad);

        char client_message[8];
        snprintf(client_message, sizeof(client_message), "P %d %d", azimuthInt, elevationInt);

        printf("%s",client_message);

        // Send the message to server:
        if(send(socket_desc, client_message, strlen(client_message), 0) < 0){
            printf("Unable to send message\n");
            return -1;
        }


        sleep(1);

    }

    close(socket_desc);
    // End CJ Code

    return 0;
}

double solveKeplersEquation(double meanAnomalyRadians, double eccentricity){
    double E_0;
    if (meanAnomalyRadians < M_PI){
        E_0 = meanAnomalyRadians + eccentricity/2;
    }
    else {
        E_0 = meanAnomalyRadians - eccentricity/2;
    }
    double E_i = E_0 - ((E_0 - eccentricity* sin(E_0) - meanAnomalyRadians)/(1-eccentricity* cos(E_0)));
    int iterations = 0;
    while((fabs(E_i-E_0) > TOLERANCE) & (iterations < MAX_ITERATIONS)){
        iterations += 1;
        E_0 = E_i;
        E_i = E_0 - ((E_0 - eccentricity* sin(E_0) - meanAnomalyRadians)/(1-eccentricity* cos(E_0)));
    }
    return E_i;
}

double calculateTrueAnomaly(double eccentricAnomaly, double eccentricity){
    return fmod(2* atan(sqrt( (1+eccentricity)/(1-eccentricity)) * tan(eccentricAnomaly/2)) + 2*M_PI,2*M_PI);
}

