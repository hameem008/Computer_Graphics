#include "2005055_classes.h"
#include "bitmap_image.hpp"

// Camera position and orientation
GLfloat eyex = 120, eyey = 120, eyez = 120;    // Camera position coordinates
GLfloat centerx = 0, centery = 0, centerz = 0; // Look-at point coordinates
GLfloat upx = 0, upy = 0, upz = 1;             // Up vector coordinates

int animation_speed = 10; // for calling timer function to simulate animation, calling after every 10 ms

// using ortho-normal vector, these three vector are perpendicular with each other
// initializing all the ortho-normal vector
Vector view_vector(centerx - eyex, centery - eyey, centerz - eyez);
Vector up_vector(0, 0, 1);
Vector right_vector(view_vector.cross(up_vector).normalize());

// image count
int image_count = 1;
// floor
Floor *flr;

// Function Declarations
void loadData();
void drawAxes();
void capture();
void initGL();
void display();
void reshapeListener(GLsizei width, GLsizei height);
void keyboardListener(unsigned char key, int x, int y);
void specialKeyListener(int key, int x, int y);

// Function to read scene.txt and populate objects and lights
void loadData()
{
    ifstream file("scene.txt");
    if (!file.is_open())
    {
        cerr << "Error: Could not open scene.txt" << endl;
        exit(1);
    }

    // Read recursion level and image resolution
    file >> recursion_level >> image_resolution;

    // Read number of objects
    int num_objects;
    file >> num_objects;

    // Read objects
    for (int i = 0; i < num_objects; ++i)
    {
        string type;
        file >> type;

        if (type == "sphere")
        {
            lld cx, cy, cz, radius;
            lld r, g, b;
            lld amb, diff, spec, refl;
            int shine;
            file >> cx >> cy >> cz >> radius;
            file >> r >> g >> b;
            file >> amb >> diff >> spec >> refl >> shine;

            Point center(cx, cy, cz);
            Color color(r, g, b);
            Coefficients coef(amb, diff, spec, refl);
            objects.push_back(new Sphere(radius, center, color, coef, shine));
        }
        else if (type == "triangle")
        {
            lld x1, y1, z1, x2, y2, z2, x3, y3, z3;
            lld r, g, b;
            lld amb, diff, spec, refl;
            int shine;
            file >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
            file >> r >> g >> b;
            file >> amb >> diff >> spec >> refl >> shine;

            Point p1(x1, y1, z1), p2(x2, y2, z2), p3(x3, y3, z3);
            Color color(r, g, b);
            Coefficients coef(amb, diff, spec, refl);
            objects.push_back(new Triangle(p1, p2, p3, color, coef, shine));
        }
        else if (type == "general")
        {
            lld A, B, C, D, E, F, G, H, I, J;
            lld bx, by, bz, bl, bw, bh;
            lld r, g, b;
            lld amb, diff, spec, refl;
            int shine;
            file >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;
            file >> bx >> by >> bz >> bl >> bw >> bh;
            file >> r >> g >> b;
            file >> amb >> diff >> spec >> refl >> shine;

            Point box_ref(bx, by, bz);
            Color color(r, g, b);
            Coefficients coef(amb, diff, spec, refl);
            objects.push_back(new GeneralQuadric(A, B, C, D, E, F, G, H, I, J, box_ref, bl, bw, bh, color, coef, shine));
        }
    }

    // Read number of point lights
    int num_point_lights;
    file >> num_point_lights;

    // Read point lights
    for (int i = 0; i < num_point_lights; ++i)
    {
        lld px, py, pz, r, g, b;
        file >> px >> py >> pz >> r >> g >> b;
        Point pos(px, py, pz);
        Color color(r, g, b);
        point_lights.push_back(PointLight(pos, color));
    }

    // Read number of spotlights
    int num_spot_lights;
    file >> num_spot_lights;

    // Read spotlights
    for (int i = 0; i < num_spot_lights; ++i)
    {
        lld px, py, pz, r, g, b, dx, dy, dz, angle;
        file >> px >> py >> pz >> r >> g >> b >> dx >> dy >> dz >> angle;
        Point pos(px, py, pz);
        Color color(r, g, b);
        Vector dir(dx, dy, dz);
        spot_lights.push_back(SpotLight(pos, color, dir, angle));
    }

    // Add floor
    Coefficients floor_coef(0.4, 0.2, 0.1, 0.3); // Example coefficients
    flr = new Floor(1000, 20, floor_coef);
    flr->loadTexture("tex_1.jpg");
    objects.push_back(flr);

    file.close();
}

void drawAxes()
{
    glLineWidth(3); // Set line thickness

    glBegin(GL_LINES);

    // X axis (red)
    glColor3f(100, 0, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(100, 0, 0);

    // Y axis (green)
    glColor3f(0, 100, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 100, 0);

    // Z axis (blue)
    glColor3f(0, 0, 100);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 100);

    glEnd();
}

// Fixed capture function
void capture()
{
    bitmap_image image(image_resolution, image_resolution);
    image.set_all_channels(0, 0, 0); // Black background

    // Calculate plane distance using correct field of view (45 degrees)
    lld viewAngle = 45.0 * M_PI / 180.0; // Convert to radians
    lld planeDistance = (image_resolution / 2.0) / tan(viewAngle / 2.0);

    Point eye(eyex, eyey, eyez);
    Vector l = view_vector.normalize();
    Vector r = right_vector.normalize();
    Vector u = up_vector.normalize();

    // Calculate topleft corner of the image plane
    Point topleft_point;
    topleft_point.x = eye.x + l.x * planeDistance - r.x * (image_resolution / 2.0) + u.x * (image_resolution / 2.0);
    topleft_point.y = eye.y + l.y * planeDistance - r.y * (image_resolution / 2.0) + u.y * (image_resolution / 2.0);
    topleft_point.z = eye.z + l.z * planeDistance - r.z * (image_resolution / 2.0) + u.z * (image_resolution / 2.0);

    // Adjust topleft to center of first pixel
    topleft_point.x += r.x * 0.5 - u.x * 0.5;
    topleft_point.y += r.y * 0.5 - u.y * 0.5;
    topleft_point.z += r.z * 0.5 - u.z * 0.5;

    // Ray tracing loop
    for (int i = 0; i < image_resolution; ++i)
    {
        for (int j = 0; j < image_resolution; ++j)
        {
            // Calculate current pixel position
            Point curPixel;
            curPixel.x = topleft_point.x + r.x * j - u.x * i;
            curPixel.y = topleft_point.y + r.y * j - u.y * i;
            curPixel.z = topleft_point.z + r.z * j - u.z * i;

            // Create ray from eye to current pixel
            Vector rayDir(curPixel.x - eye.x, curPixel.y - eye.y, curPixel.z - eye.z);
            rayDir = rayDir.normalize();
            Ray ray(eye, rayDir);

            // Find nearest intersection
            lld tMin = 1e9;
            Color pixelColor(0, 0, 0); // Default black

            // Check intersection with all objects
            for (auto *obj : objects)
            {
                Color objColor;
                lld t = obj->intersect(ray, objColor, 1);
                if (t > EPSILON && t < tMin)
                {
                    tMin = t;
                    pixelColor = objColor; // Use the object's color
                }
            }

            // Set pixel color (clamp values to [0, 255])
            int red = (int)(min((lld)255.0, max((lld)0.0, pixelColor.r * 255)));
            int green = (int)(min((lld)255.0, max((lld)0.0, pixelColor.g * 255)));
            int blue = (int)(min((lld)255.0, max((lld)0.0, pixelColor.b * 255)));

            image.set_pixel(j, i, red, green, blue);
        }
    }

    // Save image with correct naming convention
    string filename = "Output_" + to_string(image_count++) + ".bmp";
    image.save_image(filename);
    cout << "Image saved as " << filename << endl;
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

    // draw axes
    drawAxes();

    // draw objects
    for (auto *obj : objects)
    {
        if (obj != nullptr)
        {
            obj->draw();
        }
    }

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
    lld v = 0.5 * 2; // Movement increment

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
        Vector term1 = view_vector.scalar_mul(cos(theta));
        Vector term2 = up_vector.cross(view_vector).scalar_mul(sin(theta));
        Vector term3 = up_vector.scalar_mul(up_vector.dot(view_vector) * (1 - cos(theta)));
        Vector rot = term1.add(term2).add(term3);

        centerx = eyex + rot.x, centery = eyey + rot.y, centerz = eyez + rot.z;
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
        Vector term1 = view_vector.scalar_mul(cos(theta));
        Vector term2 = up_vector.cross(view_vector).scalar_mul(sin(theta));
        Vector term3 = up_vector.scalar_mul(up_vector.dot(view_vector) * (1 - cos(theta)));
        Vector rot = term1.add(term2).add(term3);

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
        Vector term1 = view_vector.scalar_mul(cos(theta));
        Vector term2 = right_vector.cross(view_vector).scalar_mul(sin(theta));
        Vector term3 = right_vector.scalar_mul(right_vector.dot(view_vector) * (1 - cos(theta)));
        Vector rot = term1.add(term2).add(term3);

        centerx = eyex + rot.x, centery = eyey + rot.y, centerz = eyez + rot.z;
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
        Vector term1 = view_vector.scalar_mul(cos(theta));
        Vector term2 = right_vector.cross(view_vector).scalar_mul(sin(theta));
        Vector term3 = right_vector.scalar_mul(right_vector.dot(view_vector) * (1 - cos(theta)));
        Vector rot = term1.add(term2).add(term3);

        centerx = eyex + rot.x, centery = eyey + rot.y, centerz = eyez + rot.z;
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
        lld theta = (v / 100) * M_PI / 180;
        Vector term1 = up_vector.scalar_mul(cos(theta));
        Vector term2 = view_vector.cross(up_vector).scalar_mul(sin(theta));
        Vector term3 = view_vector.scalar_mul(view_vector.dot((up_vector)) * (1 - cos(theta)));
        Vector rot = term1.add(term2).add(term3);

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
        lld theta = -(v / 100) * M_PI / 180;
        Vector term1 = up_vector.scalar_mul(cos(theta));
        Vector term2 = view_vector.cross(up_vector).scalar_mul(sin(theta));
        Vector term3 = view_vector.scalar_mul(view_vector.dot((up_vector)) * (1 - cos(theta)));
        Vector rot = term1.add(term2).add(term3);

        up_vector.update(rot.normalize());
        upx = up_vector.x, upy = up_vector.y, upz = up_vector.z;
        right_vector.update(view_vector.cross(up_vector).normalize());
        break;
    }

    case 'w':
    {
        // Move upward without changing reference point
        eyez += v;
        view_vector.z = centerz - eyez;
        up_vector.update(right_vector.cross(view_vector).normalize());
        upx = up_vector.x, upy = up_vector.y, upz = up_vector.z;
        break;
    }

    case 's':
    {
        // Move doqnward without changing reference point
        eyez -= v;
        view_vector.z = centerz - eyez;
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

    case '0':
        capture();
        break;
    case 't':
        flr->useTexture = !flr->useTexture;
        if (flr->useTexture)
        {
            cout << "Texture\n";
        }
        else
        {
            cout << "Chessboard\n";
        }
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
    lld v = 0.5 * 2; // Movement increment

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
        Vector view_temp(view_vector.normalize());
        eyex += view_temp.x * v, eyey += view_temp.y * v, eyez += view_temp.z * v;
        centerx += view_temp.x * v, centery += view_temp.y * v, centerz += view_temp.z * v;
        break;
    }
    case GLUT_KEY_DOWN:
    {
        // move the cammera to the backward
        // moving to the opposite direction of the view vector
        Vector view_norm(view_vector.normalize());
        eyex -= view_norm.x * v, eyey -= view_norm.y * v, eyez -= view_norm.z * v;
        centerx -= view_norm.x * v, centery -= view_norm.y * v, centerz -= view_norm.z * v;
        break;
    }
    case GLUT_KEY_PAGE_UP:
    {
        // move the cammera to the upward
        eyez += v;
        centerz += v;
        break;
    }
    case GLUT_KEY_PAGE_DOWN:
    {
        // move the cammera to the downward
        eyez -= v;
        centerz -= v;
        break;
    }
    }

    glutPostRedisplay(); // Request a screen refresh
}

void timerFunction(int value)
{

    glutPostRedisplay();

    glutTimerFunc(animation_speed, timerFunction, 0);
}

int main(int argc, char **argv)
{
    // load scene
    loadData();

    // Initialize GLUT
    glutInit(&argc, argv);

    // Configure display mode and window
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(640, 640);
    glutInitWindowPosition(50, 50);
    glutCreateWindow("Ray Tracing");

    // Register callback functions
    glutDisplayFunc(display);
    glutReshapeFunc(reshapeListener);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    // glutTimerFunc(animation_speed, timerFunction, 0);

    // Initialize OpenGL settings
    initGL();

    // Enter the GLUT event loop
    glutMainLoop();

    return 0;
}