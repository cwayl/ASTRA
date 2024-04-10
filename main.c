/*
 * Active Satellite Radio Tracking Antenna (ASTRA)
 * Coded by:
 * Cooper Gowan
 * CJ Wayland
 *
 * References to Equations are from:
 * Orbital Mechanics for Engineering Students (Forth Edition)
 *   by Howard d. Curtis
 *
 */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include "orbital_Numbers.h"
#include <stdlib.h>

#define ARRAYSIZE(a) (sizeof(a) / sizeof((a)[0]))
#define NUM_ROWS(a) ARRAYSIZE(a)
#define NUM_COLS(a) ARRAYSIZE((a)[0])
#define MAX_ITERATIONS 100 // This is the maximum number of iterations used in iterative calculations
#define TOLERANCE 1E-10 // This is the tolerance to which iterative calculations are computed
#define NUM_ROWS_Q 3 // Number of rows in Q
#define NUM_COLS_Q 3 // Number of columns in Q
#define NUM_ROWS_VECTOR 3 // Number of rows in a vector
#define NUM_COLS_VECTOR 1 // Number of columns in a vector

// Function prototypes
double solveKeplersEquation(double meanAnomaly, double eccentricity);
double calculateTrueAnomaly(double eccentricAnomaly, double eccentricity);
void string_select(const char *s, int index_start, int index_end , char *output, int size);
void moveBuffer(char *buffer, int size);

int main() {

    /*
     * The determination of where to point is broken down into ten steps.
     *
     * 1. Input TLE/Station Info
     * 2. Calculate h (specific angular momentum)
     * 3. Calculate True Anomaly
     * 4. Calculate position vector in perifocal frame
     * 5. Rotate position vector to geo equatorial frame
     *      5a: Calculate the Rotation Matrix
     *      5b: Multiply the position vector and Rotation Matrix
     * 6. Calculate sidereal time
     * 7. Calculate station position in geo equatorial frame
     * 8. Calculate relative position of satellite to station
     * 9. Calculate Azimuth and Elevation
     * 10. Check if Elevation is above horizon
     */

    // Step 1: Input TLE/Station Info

    // Define variables for user input
    double latitude, longitude, altitude;
    int epochYear;
    double epoch, inclination, raan, eccentricity, argumentOfPerigee, TLE_meanAnomaly, meanMotion;
    double prevSize = 4000;
    char buffer[100];
    char username[72];
    char password[72];
    char latitudeString[50];
    char longitudeString[50];
    char altitudeString[50];
    char bashCommand[200];
    FILE *Credentials;
    FILE *Station;
    FILE *TLE;
//    FILE *OUTPUT;
    buffer[0] = 'a';

    while (buffer[0] != '1') {
        system("bash /home/astra/ASTRA/Start-Menu.sh");
        moveBuffer(buffer, sizeof buffer);
        if (buffer[0] == '2') {
            // Setup Menu
            while (buffer[0] != '4') {
                system("bash /home/astra/ASTRA/Setup-Menu.sh");
                moveBuffer(buffer, sizeof buffer);
                if (buffer[0] == '1') {
                    // Station Coordinates

                    Station = fopen("Station.txt", "r");
                    fgets(latitudeString,sizeof latitudeString, Station);
                    fgets(longitudeString,sizeof longitudeString, Station);
                    fgets(altitudeString,sizeof altitudeString, Station);
                    fclose(Station);

                    snprintf(bashCommand, sizeof(bashCommand),
                             "bash /home/astra/ASTRA/Lat.sh %s",
                             latitudeString);
                    system(bashCommand);
                    moveBuffer(buffer,sizeof buffer);

                    Station = fopen("Station.txt", "w");
                    fprintf(Station, "%s", buffer);
                    fprintf(Station, "\n");
                    fclose(Station);

                    snprintf(bashCommand, sizeof(bashCommand),
                             "bash /home/astra/ASTRA/Long.sh %s",
                             longitudeString);
                    system(bashCommand);
                    moveBuffer(buffer,sizeof buffer);

                    Station = fopen("Station.txt", "ab");
                    fprintf(Station, "%s", buffer);
                    fprintf(Station, "\n");
                    fclose(Station);

                    snprintf(bashCommand, sizeof(bashCommand),
                             "bash /home/astra/ASTRA/Altitude.sh %s",
                             altitudeString);
                    system(bashCommand);
                    moveBuffer(buffer,sizeof buffer);

                    Station = fopen("Station.txt", "ab");
                    fprintf(Station, "%s", buffer);
                    fclose(Station);

                } else if (buffer[0] == '2') {
                    // Space-Track.com Credentials

                    Credentials = fopen("Credentials.txt", "r");
                    fgets(username, sizeof username, Credentials);
                    fgets(password, sizeof password, Credentials);
                    fclose(Credentials);

                    // Remove extra character form username variable
                    char *e;
                    int index;
                    e = strchr(username, '\n');
                    index = (int) (e - username);
                    username[index] = '\0';

                    snprintf(bashCommand, sizeof(bashCommand),
                             "bash /home/astra/ASTRA/Username.sh %s",
                             username);
                    system(bashCommand);
                    moveBuffer(buffer,sizeof buffer);

                    Credentials = fopen("Credentials.txt", "w");
                    fprintf(Credentials, "%s", buffer);
                    fprintf(Credentials, "\n");
                    fclose(Credentials);

                    snprintf(bashCommand, sizeof(bashCommand),
                             "bash /home/astra/ASTRA/Password.sh %s",
                             password);
                    system(bashCommand);
                    moveBuffer(buffer, sizeof buffer);

                    Credentials = fopen("Credentials.txt", "ab");
                    fprintf(Credentials, "%s", buffer);
                    fclose(Credentials);

                }
            }
        }
    }
    // Start ASTRA
    system("bash /home/astra/ASTRA/Manual-Select.sh");
    moveBuffer(buffer, sizeof buffer);
    // Get lat long and alt from Station file
    Station = fopen("Station.txt", "r");
    fgets(latitudeString,sizeof latitudeString, Station);
    fgets(longitudeString,sizeof longitudeString, Station);
    fgets(altitudeString,sizeof altitudeString, Station);
    fclose(Station);
    latitude = atof(latitudeString);
    longitude = atof(longitudeString);
    altitude = atof(altitudeString);

    // Convert latitude to radians for future calculations
    latitude = deg2rad(latitude);

    if (buffer[0] == '1') {
        // Use NORAD ID (Space-Track.com)

        system("bash /home/astra/ASTRA/NORAD.sh");
        moveBuffer(buffer, sizeof buffer);

        // Get the NORAD ID of the satellite
        char cookie_Command[150];
        char TLE_Command[150];

        // Get username and password form Credentials file
        Credentials = fopen("Credentials.txt", "r");
        fgets(username, sizeof username, Credentials);
        fgets(password, sizeof password, Credentials);
        fclose(Credentials);

        // Remove extra character form username variable
        char *e;
        int index;
        e = strchr(username, '\n');
        index = (int) (e - username);
        username[index] = '\0';

        // Get cookie using username and login
        snprintf(cookie_Command, sizeof cookie_Command,
                 "curl -c cookies.txt -b cookies.txt https://www.space-track.org/ajaxauth/login -d 'identity=%s&password=%s'",
                 username, password);
        system(cookie_Command);
        // Use NORAD ID to make API command
        // Puts TLE in TLE.txt
        snprintf(TLE_Command, sizeof(TLE_Command),
                 "curl --limit-rate 100K --cookie cookies.txt https://www.space-track.org/basicspacedata/query/class/gp/format/tle/NORAD_CAT_ID/%s  > TLE.txt",
                 buffer);
        system(TLE_Command);

    } else if (buffer[0] == '2'){
        // Manually Enter New TLE
        system("bash /home/astra/ASTRA/TLE1.sh");
        moveBuffer(buffer, sizeof buffer);
        TLE = fopen("TLE.txt", "w");
        fprintf(TLE, "%s", buffer);
        fprintf(TLE, "\n");
        fclose(TLE);

        system("bash /home/astra/ASTRA/TLE2.sh");
        moveBuffer(buffer, sizeof buffer);
        TLE = fopen("TLE.txt", "ab");
        fprintf(TLE, "%s", buffer);
        fclose(TLE);

    }

    // Open TLE.txt and read lines to strings
    TLE = fopen("TLE.txt", "r");
    char TLE_Line1[72];
    char TLE_Line2[72];
    fgets(TLE_Line1,sizeof TLE_Line1, TLE);
    fgets(TLE_Line2,sizeof TLE_Line2, TLE);
    fclose(TLE);

    // Copy numbers from strings to usable variables
    char splitter[13];
    string_select(TLE_Line1, 18, 19, splitter, sizeof(splitter));
    epochYear = atoi(splitter) + 2000;
    string_select(TLE_Line1, 20, 31, splitter, sizeof(splitter));
    epoch = atof(splitter);
    string_select(TLE_Line2, 8, 15, splitter, sizeof(splitter));
    inclination = atof(splitter);
    string_select(TLE_Line2, 17, 24, splitter, sizeof(splitter));
    raan = atof(splitter);
    string_select(TLE_Line2, 26, 32, splitter, sizeof(splitter));
    eccentricity = atof(splitter) * 0.0000001;
    string_select(TLE_Line2, 34, 41, splitter, sizeof(splitter));
    argumentOfPerigee = atof(splitter);
    string_select(TLE_Line2, 43, 50, splitter, sizeof(splitter));
    TLE_meanAnomaly = atof(splitter);
    string_select(TLE_Line2, 52, 62, splitter, sizeof(splitter));
    meanMotion = atof(splitter);

    // Step 2: Calculate h (specific angular momentum)

    double a; // Semi-major axis
    double h; // specific angular momentum
    double n; // Mean Motion in radians per second

    // Convert Mean Motion from revolutions per day to radians per second
    n = meanMotion * (2 * M_PI) / (24 * 3600);

    // Equation 2.83 solved for a
    // Substituting Equation 3.9 for T
    a = pow(EARTH_MU / pow(n, 2), 1.0 / 3.0);
    // Equation 2.45 solved for h
    // Substituting Equation 2.73 we can set r = r_p and theta = 0
    h = sqrt(EARTH_MU * a * (1 - pow(eccentricity,2)));

    // Perturbations
    double w_dot; // Average rate of change in argument of perigee
    double omega_dot; // Average rate of change in RAAN

    // Define the transformation matrix Q based on inclination, RAAN, and argument of perigee
    double raanRad = deg2rad(raan);
    double inclinationRad = deg2rad(inclination);
    double argumentOfPerigeeRad = deg2rad(argumentOfPerigee);

    omega_dot = -((3.0/2.0) * ((J_2 * sqrt(EARTH_MU) * pow(EARTH_RADIUS_EQUATORIAL, 2)) / (pow(a, 7.0/2) * pow((1 - pow(eccentricity, 2)), 2)))) * (cos(inclinationRad));
    w_dot = -(3.0/2.0) * ((J_2 * sqrt(EARTH_MU) * pow(EARTH_RADIUS_EQUATORIAL, 2)) /( (pow(a, 7.0/2) * pow((1 - pow(eccentricity, 2)), 2)))) * ((5.0/2.0) * pow((sin(inclinationRad)), 2) -2);

    int prevElevation, prevAzimuth;
    prevAzimuth = -1;
    prevElevation = -1;

    while(1) {

        // The current time is needed for steps three and six, so it is calculated here
        // Find UTC time
        // Structure ptr gets UTC time
        struct tm *ptr;
        time_t t;
        t = time(NULL);
        ptr = gmtime(&t);

        //Current time in hours
        double UT = (float) ptr->tm_hour + (float) ptr->tm_min / 60 + (float) ptr->tm_sec / 3600;

        // Step 3: Calculate True Anomaly

        double dayOfYear = ptr->tm_yday + 1 + UT/24;
        // Time since the TLE was observed in days
        // Assumes the TLE and current year are the same
        double delta_t = dayOfYear - epoch;
        double delta_MeanAnamoly = delta_t * meanMotion * 360;

        // Calculates current Mean anomaly
        double meanAnomaly = deg2rad(fmod(TLE_meanAnomaly + delta_MeanAnamoly, 360) );

        // Solve Kepler's Equation for Eccentric Anomaly
        double eccentricAnomaly = solveKeplersEquation(meanAnomaly, eccentricity);

        // Calculate True anomaly
        double trueAnomaly = calculateTrueAnomaly(eccentricAnomaly, eccentricity);

        // Step 4: Calculate position vector in perifocal frame
        double r_perifocal[NUM_ROWS_VECTOR][NUM_COLS_VECTOR];
        // Equation 2.119
        r_perifocal[0][0] = (h * h / EARTH_MU) * (1 / (1 + eccentricity * cos(trueAnomaly))) * cos(trueAnomaly); // x-component
        r_perifocal[1][0] = (h * h / EARTH_MU) * (1 / (1 + eccentricity * cos(trueAnomaly))) * sin(trueAnomaly); // y-component
        r_perifocal[2][0] = 0; // z-component

        //Step 5: Rotate position vector to geo equatorial frame

        double currentArgumentOfPerigeeRad, currentRaanRad;
        //TESTING
        currentArgumentOfPerigeeRad = argumentOfPerigeeRad + w_dot * delta_t * 24 * 60 * 60;
        currentRaanRad = raanRad + omega_dot * delta_t * 24 * 60 * 60;

        // Step 5a: Calculate the Rotation Matrix

        // Equation 4.49
        double R_1[NUM_ROWS_Q][NUM_COLS_Q] = {
                {cos(currentArgumentOfPerigeeRad), sin(currentArgumentOfPerigeeRad), 0},
                {-sin(currentArgumentOfPerigeeRad), cos(currentArgumentOfPerigeeRad), 0},
                {0, 0, 1},
        };
        // Equation 4.32
        double R_2[NUM_ROWS_Q][NUM_COLS_Q] = {
                {1, 0 ,0},
                {0, cos(inclinationRad), sin(inclinationRad)},
                {0, -sin(inclinationRad), cos(inclinationRad)}

        };
        // Equation 4.34
        double R_3[NUM_ROWS_Q][NUM_COLS_Q] = {
                {cos(currentRaanRad), sin(currentRaanRad), 0},
                {-sin(currentRaanRad), cos(currentRaanRad), 0},
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

        // Equation 4.49 for next four loops
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

        //      5b: Multiply the position vector and Rotation Matrix

        // Calculate the geocentric equatorial frame position vector
        // Equation 4.51
        double r_geocentric[NUM_ROWS_Q][NUM_COLS_VECTOR] = {
                {Q_Matrix_New[0][0] * r_perifocal[0][0] + Q_Matrix_New[0][1] * r_perifocal[1][0] +Q_Matrix_New[0][2] * r_perifocal[2][0]},
                {Q_Matrix_New[1][0] * r_perifocal[0][0] + Q_Matrix_New[1][1] * r_perifocal[1][0] +Q_Matrix_New[1][2] * r_perifocal[2][0]},
                {Q_Matrix_New[2][0] * r_perifocal[0][0] + Q_Matrix_New[2][1] * r_perifocal[1][0] +Q_Matrix_New[2][2] * r_perifocal[2][0]}
        };

        //Step 6: Calculate sidereal time

        // Convert UTC time into Julian time
        // Equation 5.48
        double J_O = 367 * (ptr->tm_year + 1900) - floor((7 * ((ptr->tm_year + 1900)+floor(((ptr->tm_mon + 1) + 9)/12))) / 4)
                + floor((275 * (ptr->tm_mon + 1)) / 9) + ptr->tm_mday + 1721013.5;
        // Equation 5.49
        double T_O = (J_O - 2451545) / 36525;

        // Convert Julian time into Greenwich sidereal time
        // Equation 5.50
        double greenwich_sidereal_time_0 =
                100.4606184 + 36000.77004 * T_O + 0.000387933 * pow(T_O, 2) - 2.583 * pow(10, -8) * pow(T_O, 3);
        //Equation 5.51
        double greenwich_sidereal_time = greenwich_sidereal_time_0 + 360.98564724 * UT / 24;

        // Add longitude to get local sidereal time
        // Equation 5.52
        double sidereal_time = greenwich_sidereal_time + longitude;
        sidereal_time = fmod(sidereal_time, 360);

        // Step 7: Calculate station position in geo equatorial frame

        // Convert sidereal time to radians
        sidereal_time = deg2rad(sidereal_time);

        // Calculate observer position
        double R[3][1];
        // Equation 5.56
        R[0][0] = (EARTH_RADIUS_EQUATORIAL / (sqrt(1 - (2 * EARTH_FLATTENING - pow(EARTH_FLATTENING, 2)) * pow(sin(latitude), 2))) + altitude) * cos(latitude) *
                  cos(sidereal_time);
        R[1][0] = (EARTH_RADIUS_EQUATORIAL / (sqrt(1 - (2 * EARTH_FLATTENING - pow(EARTH_FLATTENING, 2)) * pow(sin(latitude), 2))) + altitude) * cos(latitude) *
                  sin(sidereal_time);
        R[2][0] = ((EARTH_RADIUS_EQUATORIAL * pow((1 - EARTH_FLATTENING), 2)) / (sqrt(1 - (2 * EARTH_FLATTENING - pow(EARTH_FLATTENING, 2)) * pow(sin(latitude), 2))) + altitude) *
                  sin(latitude);

        // Step 8: Calculate relative position of satellite to station

        double rho_G[NUM_ROWS_VECTOR][NUM_COLS_VECTOR];
        rho_G[0][0] = r_geocentric[0][0] - R[0][0];
        rho_G[1][0] = r_geocentric[1][0] - R[1][0];
        rho_G[2][0] = r_geocentric[2][0] - R[2][0];

        // Step 9: Calculate Azimuth and Elevation

        // Make the rotation matrix
        double Q_2[NUM_ROWS_Q][NUM_COLS_Q];
        // Equation 5.62a
        Q_2[0][0] = -sin(sidereal_time);
        Q_2[0][1] = cos(sidereal_time);
        Q_2[0][2] = 0;
        Q_2[1][0] = -sin(latitude) * cos(sidereal_time);
        Q_2[1][1] = -sin(latitude) * sin(sidereal_time);
        Q_2[1][2] = cos(latitude);
        Q_2[2][0] = cos(latitude) * cos(sidereal_time);
        Q_2[2][1] = cos(latitude) * sin(sidereal_time);
        Q_2[2][2] = sin(latitude);

        double rho_R[NUM_ROWS(Q_2)][NUM_COLS(rho_G)];

        // Initializing elements of rho_R to 0.
        for (int i = 0; i < NUM_ROWS(Q_2); ++i) {
            for (int j = 0; j < NUM_COLS(rho_G); ++j) {
                rho_R[i][j] = 0;
            }
        }

        // Multiplying arrays Q_2 and rho_G to get rho_R
        for (int i = 0; i < NUM_ROWS(Q_2); ++i) {
            for (int j = 0; j < NUM_COLS(rho_G); ++j) {
                for (int k = 0; k < NUM_COLS(Q_2); ++k) {
                    rho_R[i][j] += Q_2[i][k] * rho_G[k][j];
                }
            }
        }

        // normalize rho_R
        double rho_size = sqrt(pow(rho_R[0][0], 2) + pow(rho_R[1][0], 2) + pow(rho_R[2][0], 2));
        rho_R[0][0] = rho_R[0][0] / rho_size;
        rho_R[1][0] = rho_R[1][0] / rho_size;
        rho_R[2][0] = rho_R[2][0] / rho_size;

        double rangeRate = rho_size - prevSize ;
        prevSize = rho_size;

        // convert unit vector to Azimuth and Elevation
        // Variations of Equation 5.58
        double Elevation = asin(rho_R[2][0]);
        double Azimuth = acos(rho_R[1][0] / cos(Elevation));
        if ((double) rho_R[0][0] / cos(Elevation) < 0) {
            Azimuth = 2 * M_PI - Azimuth;
        }

        // Step 10: Check if Elevation is above horizon
        char scriptMessage[40];
        int azimuthInt, elevationInt;
//      if the satellite is below the horizon move to home otherwise move to vector position
//      the current position of the rotator is checked to make sure we don't repeat commands
        if (rho_R[2][0] < 0) {
            if((prevAzimuth != 180) || (prevElevation != 0)){
                system("/home/astra/rotatorScript.sh 180 0");
                prevAzimuth = 180;
                prevElevation = 0;
            }
        } else {
            azimuthInt = (int) floor(rad2deg(Azimuth));
            elevationInt = (int) floor(rad2deg(Elevation));
            if ((prevAzimuth != azimuthInt) || (prevElevation != elevationInt)) {
                snprintf(scriptMessage, sizeof(scriptMessage), "/home/astra/rotatorScript.sh %d %d",
                         azimuthInt, elevationInt);
                system(scriptMessage);
                prevElevation = elevationInt;
                prevAzimuth = azimuthInt;
            }
        }

        char outputCommand[200];
        snprintf(outputCommand, sizeof outputCommand,
                 "bash /home/astra/ASTRA/Monitor.sh %f %f %f %f",
                 rad2deg(Azimuth),
                 rad2deg(Elevation),
                 rho_size,
                 rangeRate);
        system(outputCommand);

/*
        char outputingThings[200];
        snprintf(outputingThings, sizeof outputingThings,
                 "%f, %f, %s",
                 rad2deg(Azimuth),
                 rad2deg(Elevation),
                 ctime(&t));
        OUTPUT = fopen("testSomething.txt", "ab");
        fprintf(OUTPUT, "%s", outputingThings);
        fclose(OUTPUT);
*/

        sleep(1);
    }
    return 0;
}

double solveKeplersEquation(double meanAnomalyRadians, double eccentricity){
    double E_0;
    // Initial values from Matt Harris lecture notes page 5.12
    if (meanAnomalyRadians < M_PI){
        E_0 = meanAnomalyRadians + eccentricity/2;
    }
    else {
        E_0 = meanAnomalyRadians - eccentricity/2;
    }
    // Equation 3.17
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
    // Equation 3.13a solved for true anomaly
    return fmod(2* atan(sqrt( (1+eccentricity)/(1-eccentricity)) * tan(eccentricAnomaly/2)) + 2*M_PI,2*M_PI);
}

void string_select(const char *s, int index_start, int index_end , char *output, int size)
{
    // sets the output string to zero then copies the specified part of the string to output
    for (int i = 0; i < size; i++) {
        output[i] = 0;
    }

    for (int i = index_start; i <= index_end; i++){
        output[i - index_start] = s[i];
    }
}

void moveBuffer(char *buffer, int size){
    FILE *BUFFER;
    BUFFER = fopen("/home/astra/ASTRA/text.txt", "r");
    fgets(buffer, size, BUFFER);
    fclose(BUFFER);
}
