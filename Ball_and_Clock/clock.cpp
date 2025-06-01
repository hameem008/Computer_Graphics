#include <GL/glut.h>
#include <bits/stdc++.h>
#include <ctime>
#include <chrono>
using namespace std;
using namespace std::chrono;
#define endl '\n'

int animation_speed = 100;

void init()
{
    glClearColor(0.0f, 0.0f, 0.2f, 1.0f);
    glEnable(GL_DEPTH_TEST);
}

void drawMarker(float center_x, float center_y)
{
    float radius = 0.04;      // radius of the marker circle
    int number_of_side = 100; // drawing a polygon with 100 side, as it will look like a circle
    glBegin(GL_TRIANGLE_FAN); // for filling with color
    for (int i = 0; i < number_of_side; i++)
    {
        float angle = 2.0f * M_PI * i / number_of_side;
        float x = radius * cos(angle);
        float y = radius * sin(angle);
        glVertex3f(center_x + x, center_y + y, 0.0f);
    }
    glEnd();
}

void drawClockDial(float radius, int number_of_side)
{
    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINE_LOOP); // drawing a ploygon, as it will look like a circle
    for (int i = 0; i < number_of_side; i++)
    {
        float angle = 2.0f * M_PI * i / number_of_side;
        float x = radius * cos(angle);
        float y = radius * sin(angle);
        glVertex3f(x, y, 0.0f);
    }
    glEnd();
}

void drawClockMarks(float length, int number_of_line)
{
    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINES); // drawing 60 lines and dividing them in m1:m2 ratio
    for (int i = 0; i < number_of_line; i++)
    {
        float angle = 2.0f * M_PI * i / number_of_line;
        float x = length * cos(angle);
        float y = length * sin(angle);
        int m1 = 10, m2 = 1;
        if (i % 5 == 0) // the ratio is different for some lines
            m1 = 5;
        glVertex3f((m1 * x + m2 * 0) / (m1 + m2), (m1 * y + m2 * 0) / (m1 + m2), 0.0f);
        glVertex3f(x, y, 0.0f);
    }
    glEnd();
}

void drawClockHand(float length)
{
    // getting the system time
    auto now = system_clock::now();
    time_t t = system_clock::to_time_t(now);
    tm *localTime = localtime(&t);
    int h = (localTime->tm_hour) % 12;
    int m = localTime->tm_min;
    int s = localTime->tm_sec;
    auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;

    // drawing hour hand
    glColor3f(1.0f, 1.0f, 1.0f);
    glLineWidth(6.0f);
    glBegin(GL_LINES);
    float hour_hand_cof = 0.45;                                                // for adjustuing the length of the hour hand
    float hour_angle = 2.0f * M_PI * (h * 3600 + m * 60 + s) / (3600 * 12.0f); // hour:minute:second = 3600:60:1
    float hour_x = length * sin(hour_angle);
    float hour_y = length * cos(hour_angle);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(hour_x * hour_hand_cof, hour_y * hour_hand_cof, 0.0f);
    glEnd();
    glLineWidth(1.0f);
    // drawing hour hand marker
    glColor3f(1.0f, 1.0f, 1.0f);
    drawMarker(hour_x, hour_y);

    // minute hand
    glColor3f(1.0f, 1.0f, 1.0f);
    glLineWidth(3.0f);
    glBegin(GL_LINES);
    float minute_hand_cof = 0.7;                                    // for adjustuing the length of the minute hand
    float minute_angle = 2.0f * M_PI * (m * 60 + s) / (60 * 60.0f); // minute:second = 60:1
    float minute_x = length * sin(minute_angle);
    float minute_y = length * cos(minute_angle);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(minute_x * minute_hand_cof, minute_y * minute_hand_cof, 0.0f);
    glEnd();
    glLineWidth(1.0f);
    // drawing minute hand marker
    glColor3f(1.0f, 1.0f, 1.0f);
    drawMarker(minute_x, minute_y);

    // second hand
    glColor3f(1.0f, 0.0f, 0.0f);
    glBegin(GL_LINES);
    float second_hand_cof = 0.7;                                                      // for adjustuing the length of the second hand
    float second_angle = 2.0f * M_PI * (s * 10 + ms.count() / 100.0f) / (10 * 60.0f); // s:(ms/100) = 10:1
    float second_x = length * sin(second_angle);
    float second_y = length * cos(second_angle);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(second_x * second_hand_cof, second_y * second_hand_cof, 0.0f);
    glEnd();
    // drawing second hand marker
    glColor3f(1.0f, 0.0f, 0.0f);
    drawMarker(second_x, second_y);
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();

    gluLookAt(
        0.0, 0.0, 3.0,
        0.0, 0.0, 0.0,
        0.0, 1.0, 0.0);

    drawClockDial(1.05, 10000);
    drawClockMarks(1, 60);
    drawClockHand(1.05);

    glutSwapBuffers();
}

void reshape(int width, int height)
{
    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, (float)width / height, 1.0, 100.0);

    glMatrixMode(GL_MODELVIEW);
}

void timerFunction(int value)
{
    glutPostRedisplay();

    glutTimerFunc(animation_speed, timerFunction, 0);
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);

    glutInitWindowSize(1000, 800);
    glutInitWindowPosition(100, 100);

    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("Clock");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutTimerFunc(animation_speed, timerFunction, 0);
    init();

    cout << "Main loop" << endl;
    glutMainLoop();

    return 0;
}