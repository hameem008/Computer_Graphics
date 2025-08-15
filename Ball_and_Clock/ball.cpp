#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;
typedef long double lld;
typedef long long ll;
#include <GL/glut.h>

class vect
{
public:
    lld x, y, z;
    vect(lld a, lld b, lld c)
    {
        x = a, y = b, z = c;
    }
    vect(const vect &v)
    {
        x = v.x, y = v.y, z = v.z;
    }
    lld dot(vect v)
    {
        return (x * v.x + y * v.y + z * v.z);
    }
    vect cross(vect v)
    {
        vect ret(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
        return ret;
    }
    vect add(vect v)
    {
        vect ret(x + v.x, y + v.y, z + v.z);
        return ret;
    }
    vect scalar_mul(lld m)
    {
        vect ret(x * m, y * m, z * m);
        return ret;
    }
    lld length()
    {
        return sqrt(x * x + y * y + z * z);
    }
    vect normalize()
    {
        lld l = length();
        vect ret(x / l, y / l, z / l);
        return ret;
    }
    void update(vect v)
    {
        x = v.x, y = v.y, z = v.z;
    }
    void update(lld a, lld b, lld c)
    {
        x = a, y = b, z = c;
    }
};

// Camera position and orientation
// GLfloat eyex = 30, eyey = 0, eyez = 30;        // Camera position coordinates
// GLfloat centerx = 0, centery = 0, centerz = 0; // Look-at point coordinates
// GLfloat upx = 0, upy = 1, upz = 0;             // Up vector coordinates

// for better view
GLfloat eyex = 37.7217, eyey = 21.9431, eyez = 37.7217;           // Camera position coordinates
GLfloat centerx = 5.40567, centery = -1.31225, centerz = 5.40567; // Look-at point coordinates
GLfloat upx = -0.320683, upy = 0.891249, upz = -0.320683;         // Up vector coordinates

const lld BOX_WIDTH = 24.0f;
const lld BOX_HEIGHT = 12.0f;
lld ball_radius = 0.5f;
// the range of x coordinate value for the ball, also for collisiion ditection on yz plane
lld ball_xmin = ball_radius, ball_xmax = BOX_WIDTH - ball_radius;
// the range of y coordinate value for the ball, also for collisiion ditection on zx plane
lld ball_ymin = ball_radius, ball_ymax = BOX_HEIGHT - ball_radius;
// the range of z coordinate value for the ball, also for collisiion ditection on xy plane
lld ball_zmin = ball_radius, ball_zmax = BOX_WIDTH - ball_radius;
lld ball_x = 1.0, ball_y = ball_radius, ball_z = 1.0; // coordiante of the ball, here initialized
// for drawing velocity arrow
lld arrow_x1 = ball_x, arrow_y1 = ball_y, arrow_z1 = ball_z;
lld arrow_x2 = arrow_x1, arrow_y2 = 2 * arrow_y1, arrow_z2 = arrow_z1;
// this angle is for ball spinning
lld ball_angle = 0;

int animation_speed = 10;                  // for calling timer function to simulate animation, calling after every 10 ms
bool is_running = false;                   // is the ball is moving or paused
bool is_arrow = true;                      // for toggleable arrow
lld damping = 0.8;                         // damping coefficient for collision
lld velocity_low = 20, velocity_high = 30; // velocity range for randomization

// using ortho-normal vector, these three vector are perpendicular with each other
// initializing all the ortho-normal vector
vect view_vector(centerx - eyex, centery - eyey, centerz - eyez);
vect up_vector(0, 1, 0);
vect right_vector(view_vector.cross(up_vector).normalize());
// velocity vector for the ball
vect ball_velocity((velocity_low + velocity_high) / 2, (velocity_low + velocity_high) / 2, (velocity_low + velocity_high) / 2);
vect plane_normal(0, 1, 0); // for calculating the spinning axis for the ball
vect spin_axis(0, 0, 0);    // spnning axis of the ball

// Function Declarations
void initGL();
void display();
void reshapeListener(GLsizei width, GLsizei height);
void keyboardListener(unsigned char key, int x, int y);
void specialKeyListener(int key, int x, int y);
void drawArrowHead();
void drawArrowShaft();
void drawBox();

// for custom randomization in a given range
ll random(lld low, lld high)
{
    ll lowll = static_cast<ll>(ceil(low));
    ll highll = static_cast<ll>(floor(high));
    return lowll + rand() % (highll - lowll + 1);
}

/**
 * Initialize OpenGL settings
 * Sets up background color and enables depth testing
 */
void initGL()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black background
    glEnable(GL_DEPTH_TEST);              // Enable depth testing for z-culling
}

void display()
{
    // Clear color and depth buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Set up the model-view matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Position camera using the eye, center and up vectors
    gluLookAt(eyex, eyey, eyez,          // Camera position
              centerx, centery, centerz, // Look-at point
              upx, upy, upz);            // Up vector

    // Draw objects based on visibility flags
    drawBox();

    // Swap buffers (double buffering)
    glutSwapBuffers();
}

void reshapeListener(GLsizei width, GLsizei height)
{
    // Prevent division by zero
    if (height == 0)
        height = 1;

    // Calculate aspect ratio
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set viewport to cover entire window
    glViewport(0, 0, width, height);

    // Set up perspective projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // 45-degree field of view, aspect ratio, near and far clipping planes
    gluPerspective(45.0f, aspect, 0.1f, 10000.0f);
    // changed the far clipping plane
}

void keyboardListener(unsigned char key, int x, int y)
{
    lld v = 0.8; // Movement increment

    // rodrigues formulla
    // a vector v will be rotated by an angle θ aroud axis k
    // vr = vcos(θ) + (k x v)sin(θ) + k(k . v)(1 - cos(θ))

    switch (key)
    {
    // Camera Position Controls
    case '1':
    {
        // look left (yaw)
        // rodrigues formulla is applied here
        // rotation axis up_vector
        // rotation vector view_vector
        lld theta = v * M_PI / 180;
        vect term1 = view_vector.scalar_mul(cos(theta));
        vect term2 = up_vector.cross(view_vector).scalar_mul(sin(theta));
        vect term3 = up_vector.scalar_mul(up_vector.dot(view_vector) * (1 - cos(theta)));
        vect rot = term1.add(term2).add(term3);

        centerx = eyex + rot.x;
        centery = eyey + rot.y;
        centerz = eyez + rot.z;
        view_vector.update(rot);
        right_vector.update(view_vector.cross(up_vector).normalize());
        break;
    }

    case '2':
    {
        // look right (yaw)
        // rodrigues formulla is applied here
        // rotation axis up_vector
        // rotation vector view_vector
        lld theta = -v * M_PI / 180;
        vect term1 = view_vector.scalar_mul(cos(theta));
        vect term2 = up_vector.cross(view_vector).scalar_mul(sin(theta));
        vect term3 = up_vector.scalar_mul(up_vector.dot(view_vector) * (1 - cos(theta)));
        vect rot = term1.add(term2).add(term3);

        centerx = eyex + rot.x, centery = eyey + rot.y, centerz = eyez + rot.z;
        view_vector.update(rot);
        right_vector.update(view_vector.cross(up_vector).normalize());
        break;
    }

    case '3':
    {
        // Look up (pitch)
        // rodrigues formulla is applied here
        // rotation axis right_vector
        // rotation vector view_vector
        lld theta = v * M_PI / 180;
        vect term1 = view_vector.scalar_mul(cos(theta));
        vect term2 = right_vector.cross(view_vector).scalar_mul(sin(theta));
        vect term3 = right_vector.scalar_mul(right_vector.dot(view_vector) * (1 - cos(theta)));
        vect rot = term1.add(term2).add(term3);

        centerx = eyex + rot.x;
        centery = eyey + rot.y;
        centerz = eyez + rot.z;
        view_vector.update(rot);
        up_vector.update(right_vector.cross(view_vector).normalize());
        upx = up_vector.x, upy = up_vector.y, upz = up_vector.z;
        break;
    }

    case '4':
    {
        // Look down (pitch)
        // rodrigues formulla is applied here
        // rotation axis right_vector
        // rotation vector view_vector
        lld theta = -v * M_PI / 180;
        vect term1 = view_vector.scalar_mul(cos(theta));
        vect term2 = right_vector.cross(view_vector).scalar_mul(sin(theta));
        vect term3 = right_vector.scalar_mul(right_vector.dot(view_vector) * (1 - cos(theta)));
        vect rot = term1.add(term2).add(term3);

        centerx = eyex + rot.x;
        centery = eyey + rot.y;
        centerz = eyez + rot.z;
        view_vector.update(rot);
        up_vector.update(right_vector.cross(view_vector).normalize());
        upx = up_vector.x, upy = up_vector.y, upz = up_vector.z;
        break;
    }

    case '5':
    {
        // Tilt clockwise (roll)
        // rodrigues formulla is applied here
        // rotation axis view_vector
        // rotation vector up_vector
        lld theta = (v / 10) * M_PI / 180;
        vect term1 = up_vector.scalar_mul(cos(theta));
        vect term2 = view_vector.cross(up_vector).scalar_mul(sin(theta));
        vect term3 = view_vector.scalar_mul(view_vector.dot((up_vector)) * (1 - cos(theta)));
        vect rot = term1.add(term2).add(term3);

        up_vector.update(rot.normalize());
        upx = up_vector.x, upy = up_vector.y, upz = up_vector.z;
        right_vector.update(view_vector.cross(up_vector).normalize());
        break;
    }

    case '6':
    {
        // Tilt anticlockwise (roll)
        // rodrigues formulla is applied here
        // rotation axis view_vector
        // rotation vector up_vector
        lld theta = -(v / 10) * M_PI / 180;
        vect term1 = up_vector.scalar_mul(cos(theta));
        vect term2 = view_vector.cross(up_vector).scalar_mul(sin(theta));
        vect term3 = view_vector.scalar_mul(view_vector.dot((up_vector)) * (1 - cos(theta)));
        vect rot = term1.add(term2).add(term3);

        up_vector.update(rot.normalize());
        upx = up_vector.x, upy = up_vector.y, upz = up_vector.z;
        right_vector.update(view_vector.cross(up_vector).normalize());
        break;
    }

    case 'w':
    {
        // Move upward without changing reference point
        eyey += v;
        view_vector.y = centery - eyey;
        up_vector.update(right_vector.cross(view_vector).normalize());
        upx = up_vector.x, upy = up_vector.y, upz = up_vector.z;
        break;
    }

    case 's':
    {
        // Move doqnward without changing reference point
        eyey -= v;
        view_vector.y = centery - eyey;
        up_vector.update(right_vector.cross(view_vector).normalize());
        upx = up_vector.x, upy = up_vector.y, upz = up_vector.z;
        break;
    }

    // for debugging perpose
    case 'h':
        cout << "centerx: " << centerx << " " << "centery: " << centery << " " << "centerz: " << centerz << endl;
        cout << "eyex: " << eyex << " " << "eyey: " << eyey << " " << "eyez: " << eyez << endl;
        cout << "upx: " << upx << " " << "upy: " << upy << " " << "upz: " << upz << endl;
        break;

    // for reseting the position and the velocity of the ball
    case 'r':
        ball_x = random(ball_xmin, ball_xmax), ball_y = ball_ymin, ball_z = random(ball_zmin, ball_zmax);
        ball_velocity.update(random(velocity_low, velocity_high), random(velocity_low, velocity_high), random(velocity_low, velocity_high));
        arrow_x1 = ball_x, arrow_y1 = ball_y, arrow_z1 = ball_z;
        arrow_x2 = arrow_x1, arrow_y2 = 2 * arrow_y1, arrow_z2 = arrow_z1;
        break;

    // for increasing ball velocity
    case '+':
        if (!is_running)
            ball_velocity.update(ball_velocity.x / 0.8, ball_velocity.y / 0.8, ball_velocity.z / 0.8);
        break;

    // for decreasing ball velocity
    case '-':
        if (!is_running)
            ball_velocity.update(ball_velocity.x * 0.8, ball_velocity.y * 0.8, ball_velocity.z * 0.8);
        break;

    // pausing the ball
    case ' ':
        is_running = !is_running;
        break;

    // toggleable velocity arror
    case 'v':
        is_arrow = !is_arrow;
        break;

    // --- Program Control ---
    case 27:
        exit(0);
        break; // ESC key: exit program
    }

    glutPostRedisplay(); // Request a screen refresh
}

void specialKeyListener(int key, int x, int y)
{
    lld v = 0.15; // Movement increment

    switch (key)
    {
    case GLUT_KEY_LEFT:
    {
        // move the cammera to the left
        // moving to the direction of the right vector
        eyex -= right_vector.x * v, eyey -= right_vector.y * v, eyez -= right_vector.z * v;
        centerx -= right_vector.x * v, centery -= right_vector.y * v, centerz -= right_vector.z * v;
        break;
    }
    case GLUT_KEY_RIGHT:
    {
        // move the cammera to the right
        // moving to the opposite direction of the right vector
        eyex += right_vector.x * v, eyey += right_vector.y * v, eyez += right_vector.z * v;
        centerx += right_vector.x * v, centery += right_vector.y * v, centerz += right_vector.z * v;
        break;
    }
    case GLUT_KEY_UP:
    {
        // move the cammera to the forward
        // moving to the direction of the view vector
        vect view_temp(view_vector.normalize());
        eyex += view_temp.x * v, eyey += view_temp.y * v, eyez += view_temp.z * v;
        centerx += view_temp.x * v, centery += view_temp.y * v, centerz += view_temp.z * v;
        break;
    }
    case GLUT_KEY_DOWN:
    {
        // move the cammera to the backward
        // moving to the opposite direction of the view vector
        vect view_norm(view_vector.normalize());
        eyex -= view_norm.x * v, eyey -= view_norm.y * v, eyez -= view_norm.z * v;
        centerx -= view_norm.x * v, centery -= view_norm.y * v, centerz -= view_norm.z * v;
        break;
    }
    case GLUT_KEY_PAGE_UP:
    {
        // move the cammera to the upward
        eyey += v;
        centery += v;
        break;
    }
    case GLUT_KEY_PAGE_DOWN:
    {
        // move the cammera to the downward
        eyey -= v;
        centery -= v;
        break;
    }
    }

    glutPostRedisplay(); // Request a screen refresh
}

void drawArrowHead(lld cx, lld cy, lld cz, lld radius, lld height, lld px, lld py, lld pz)
{
    // Direction vector from base to peak (tip of the cone)
    lld dx = px - cx;
    lld dy = py - cy;
    lld dz = pz - cz;

    // Length of the vector (distance from base to peak)
    lld length = sqrt(dx * dx + dy * dy + dz * dz);

    // Normalize the direction vector
    dx /= length;
    dy /= length;
    dz /= length;

    // Create perpendicular vectors for drawing the circular base
    // We can take a cross-product of the direction vector with any vector not parallel to it.
    lld perp_x = -dy, perp_y = dx, perp_z = 0; // Perpendicular vector in the x-y plane
    lld perp_length = sqrt(perp_x * perp_x + perp_z * perp_z);
    if (perp_length > 0)
    {
        perp_x /= perp_length;
        perp_z /= perp_length;
    }

    // Second perpendicular vector in the z-y plane
    lld side_x = dz * perp_y - dy * perp_z;
    lld side_y = dx * perp_z - dz * perp_x;
    lld side_z = dy * perp_x - dx * perp_y;

    // Normalize the second perpendicular vector
    lld sideLength = sqrt(side_x * side_x + side_y * side_y + side_z * side_z);
    side_x /= sideLength;
    side_y /= sideLength;
    side_z /= sideLength;

    int numSegments = 20; // Number of segments for the cone base circle

    // Draw the base of the cone using a triangle fan
    glColor3f(1.0f, 0.0f, 0.0f); // Red color for the cone base
    glBegin(GL_TRIANGLE_FAN);

    // Draw the peak of the cone
    glVertex3f(px, py, pz);

    // Draw the circular base of the cone
    for (int i = 0; i <= numSegments; i++)
    {
        lld angle = 2.0f * M_PI * i / numSegments; // Angle around the circle
        // Compute the coordinates around the circular base
        lld x = cx + radius * cos(angle) * perp_x + radius * sin(angle) * side_x;
        lld y = cy + radius * cos(angle) * perp_y + radius * sin(angle) * side_y;
        lld z = cz + radius * cos(angle) * perp_z + radius * sin(angle) * side_z;

        glVertex3f(x, y, z);
    }
    glEnd();

    // Optionally, you can draw the sides of the cone (lines from base points to the peak)
    glBegin(GL_LINES);
    for (int i = 0; i <= numSegments; i++)
    {
        lld angle = 2.0f * M_PI * i / numSegments; // Angle for each base point
        // Base points (x, y, z) for the circular base
        lld x = cx + radius * cos(angle) * perp_x + radius * sin(angle) * side_x;
        lld y = cy + radius * cos(angle) * perp_y + radius * sin(angle) * side_y;
        lld z = cz + radius * cos(angle) * perp_z + radius * sin(angle) * side_z;

        glVertex3f(x, y, z);    // Base point
        glVertex3f(px, py, pz); // Peak of the cone
    }
    glEnd();
}

void drawArrowShaft(float x1, lld y1, lld z1, lld x2, lld y2, lld z2)
{
    lld arrow_head_length = 0.3f;
    lld arrow_line_length = 2.0f;
    vect arrow(x2 - x1, y2 - y1, z2 - z1);
    arrow.update(arrow.normalize());
    x2 = x1 + arrow.x * arrow_line_length;
    y2 = y1 + arrow.y * arrow_line_length;
    z2 = z1 + arrow.z * arrow_line_length;

    // Set line thickness
    glLineWidth(3.0f); // Try 2.0f to 5.0f for different thickness levels

    // Draw the shaft of the arrow
    glColor3f(1.0f, 0.0f, 0.0f); // Red arrow
    glBegin(GL_LINES);
    glVertex3f(x1, y1, z1);
    glVertex3f(x2, y2, z2);
    glEnd();

    // Reset to default line width if necessary
    glLineWidth(1.0f);

    drawArrowHead(x2 - arrow.x * arrow_head_length, y2 - arrow.y * arrow_head_length, z2 - arrow.z * arrow_head_length, 0.1, 0.5, x2, y2, z2);
}

void drawBall()
{
    if (is_arrow)
        drawArrowShaft(arrow_x1, arrow_y1, arrow_z1, arrow_x2, arrow_y2, arrow_z2);

    int slices = 200;
    int stacks = 100;
    lld radius = 0.5f;

    glPushMatrix();                       // Save the current transformation
    glTranslatef(ball_x, ball_y, ball_z); // Move the ball to the given position
    glRotatef(ball_angle, spin_axis.x, spin_axis.y, spin_axis.z);

    for (int i = 0; i < stacks; ++i)
    {
        lld theta1 = (float)i / stacks * M_PI;
        lld theta2 = (float)(i + 1) / stacks * M_PI;

        glBegin(GL_QUAD_STRIP);
        for (int j = 0; j <= slices; ++j)
        {
            lld phi = (float)j / slices * 2.0f * M_PI;

            // Compute vertex positions
            lld x1 = radius * sin(theta1) * cos(phi);
            lld y1 = radius * cos(theta1);
            lld z1 = radius * sin(theta1) * sin(phi);

            lld x2 = radius * sin(theta2) * cos(phi);
            lld y2 = radius * cos(theta2);
            lld z2 = radius * sin(theta2) * sin(phi);

            // Select color based on longitude (phi)
            int band = (j * 6) / slices; // 6 vertical bands

            if (band == 0)
                glColor3f(0.9f, 0.4f, 0.4f); // Soft Red
            else if (band == 1)
                glColor3f(0.4f, 0.6f, 0.9f); // Soft Blue
            else if (band == 2)
                glColor3f(0.4f, 0.9f, 0.5f); // Soft Green
            else if (band == 3)
                glColor3f(0.9f, 0.8f, 0.4f); // Soft Yellow
            else if (band == 4)
                glColor3f(0.8f, 0.4f, 0.9f); // Soft Purple
            else if (band == 5)
                glColor3f(0.4f, 0.9f, 0.9f); // Soft Cyan

            glNormal3f(x1 / radius, y1 / radius, z1 / radius);
            glVertex3f(x1, y1, z1);

            glNormal3f(x2 / radius, y2 / radius, z2 / radius);
            glVertex3f(x2, y2, z2);
        }
        glEnd();
    }
    glPopMatrix(); // Restore the previous transformation
}

void drawBox()
{
    lld boxWidth = BOX_WIDTH;   // width and depth (for square floor/ceiling)
    lld boxHeight = BOX_HEIGHT; // height (for rectangular walls)
    int tilesPerSide = 12;      // chessboard floor
    lld tileSize = boxWidth / tilesPerSide;

    glBegin(GL_QUADS);

    // Top face (ceiling) - Light Blue
    glColor3f(0.6f, 0.8f, 1.0f);
    glVertex3f(0.0f, boxHeight, 0.0f);
    glVertex3f(boxWidth, boxHeight, 0.0f);
    glVertex3f(boxWidth, boxHeight, boxWidth);
    glVertex3f(0.0f, boxHeight, boxWidth);

    // Bottom face (floor) - Chessboard pattern
    for (int i = 0; i < tilesPerSide; ++i)
    {
        for (int j = 0; j < tilesPerSide; ++j)
        {
            if ((i + j) % 2 == 0)
                glColor3f(1.0f, 1.0f, 1.0f); // White
            else
                glColor3f(0.0f, 0.0f, 0.0f); // Black

            lld x_start = i * tileSize;
            lld z_start = j * tileSize;

            glVertex3f(x_start, 0.0f, z_start);
            glVertex3f(x_start + tileSize, 0.0f, z_start);
            glVertex3f(x_start + tileSize, 0.0f, z_start + tileSize);
            glVertex3f(x_start, 0.0f, z_start + tileSize);
        }
    }

    // Front face (wall) - Light Coral
    glColor3f(1.0f, 0.6f, 0.6f);
    glVertex3f(0.0f, boxHeight, 0.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(boxWidth, 0.0f, 0.0f);
    glVertex3f(boxWidth, boxHeight, 0.0f);

    // Back face (wall) - Light Yellow
    glColor3f(1.0f, 1.0f, 0.6f);
    glVertex3f(boxWidth, boxHeight, boxWidth);
    glVertex3f(boxWidth, 0.0f, boxWidth);
    glVertex3f(0.0f, 0.0f, boxWidth);
    glVertex3f(0.0f, boxHeight, boxWidth);

    // Left face (wall) - Light Green
    glColor3f(0.6f, 1.0f, 0.6f);
    glVertex3f(0.0f, boxHeight, boxWidth);
    glVertex3f(0.0f, 0.0f, boxWidth);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, boxHeight, 0.0f);

    // Right face (wall) - Light Purple
    glColor3f(0.8f, 0.6f, 1.0f);
    glVertex3f(boxWidth, boxHeight, 0.0f);
    glVertex3f(boxWidth, 0.0f, 0.0f);
    glVertex3f(boxWidth, 0.0f, boxWidth);
    glVertex3f(boxWidth, boxHeight, boxWidth);

    glEnd();

    drawBall();
}

void timerFunction(int value)
{
    if (is_running)
    {
        lld dt = 0.01; // as the function is called after every 10 ms
        // storing the previous ball position to draw arrow
        arrow_x1 = ball_x, arrow_y1 = ball_y, arrow_z1 = ball_z;
        // applying the velocity to the ball position
        ball_x = ball_x + ball_velocity.x * dt, ball_y = ball_y + ball_velocity.y * dt, ball_z = ball_z + ball_velocity.z * dt;
        // storing the current ball position to draw arrow
        arrow_x2 = ball_x, arrow_y2 = ball_y, arrow_z2 = ball_z;
        // collision condition
        if (ball_x < ball_xmin || ball_x > ball_xmax)
        {
            if (ball_x < ball_xmin)
                ball_x = ball_xmin;
            else if (ball_x > ball_xmax)
                ball_x = ball_xmax;
            plane_normal.update(1, 0, 0); // changing the plane normal for side spin
            ball_velocity.x = -damping * ball_velocity.x;
        }
        if (ball_y <= ball_ymin || ball_y >= ball_ymax)
        {
            if (ball_y < ball_ymin)
                ball_y = ball_ymin;
            else if (ball_y > ball_ymax)
                ball_y = ball_ymax;
            plane_normal.update(0, 1, 0); // // changing the plane normal for forward spin
            ball_velocity.y = -damping * ball_velocity.y;
        }
        if (ball_z <= ball_zmin || ball_z >= ball_zmax)
        {
            if (ball_z < ball_zmin)
                ball_z = ball_zmin;
            else if (ball_z > ball_zmax)
                ball_z = ball_zmax;
            plane_normal.update(0, 0, 1); // changing the plane normal for side spin
            ball_velocity.z = -damping * ball_velocity.z;
        }
        ball_velocity.y -= 9.8 * dt;
        // calculating rolling angle based on displacement
        ball_angle += 360 * dt * ball_velocity.length() / (2 * M_PI * ball_radius);
        spin_axis.update(plane_normal.cross(ball_velocity)); // updating spinning axis vector as velocity and point normal is variable
    }

    glutPostRedisplay();

    glutTimerFunc(animation_speed, timerFunction, 0);
}

int main(int argc, char **argv)
{
    // Initialize GLUT
    glutInit(&argc, argv);

    // Configure display mode and window
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(640, 640);
    glutInitWindowPosition(50, 50);
    glutCreateWindow("Ball");

    // Register callback functions
    glutDisplayFunc(display);
    glutReshapeFunc(reshapeListener);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutTimerFunc(animation_speed, timerFunction, 0);

    // Initialize OpenGL settings
    initGL();

    // Enter the GLUT event loop
    glutMainLoop();

    return 0;
}
