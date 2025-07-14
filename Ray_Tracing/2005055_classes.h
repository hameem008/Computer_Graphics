#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#pragma GCC optimize("O3")
#pragma GCC optimize("fast-math")
#pragma GCC optimize("unroll-loops")
#include <bits/stdc++.h>
#include <GL/glut.h>
using namespace std;
typedef double lld;

const lld EPSILON = 1e-6;

class Point
{
public:
    lld x, y, z;
    Point() { x = y = z = 0; }
    Point(lld x, lld y, lld z) { this->x = x, this->y = y, this->z = z; }
    void update(Point &p) { this->x = p.x, this->y = p.y, this->z = p.z; }
    void update(lld x, lld y, lld z) { this->x = x, this->y = y, this->z = z; }
};

class Vector
{
public:
    lld x, y, z;
    Vector() { x = y = z = 0; }
    Vector(lld x, lld y, lld z) { this->x = x, this->y = y, this->z = z; }
    Vector(const Vector &v) { this->x = v.x, this->y = v.y, this->z = v.z; }
    lld dot(Vector &v) { return (x * v.x + y * v.y + z * v.z); }
    Vector cross(Vector &v) { return Vector(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }
    Vector add(Vector &v) { return Vector(x + v.x, y + v.y, z + v.z); }
    Vector sub(Vector &v) { return Vector(x - v.x, y - v.y, z - v.z); }
    Vector scalar_mul(lld m) { return Vector(x * m, y * m, z * m); }
    lld length() { return sqrt(x * x + y * y + z * z); }
    Vector normalize()
    {
        lld l = length();
        return Vector(x / l, y / l, z / l);
    }
    void update(Vector v) { this->x = v.x, this->y = v.y, this->z = v.z; }
    void update(lld x, lld y, lld z) { this->x = x, this->y = y, this->z = z; }
};

class Color
{
public:
    lld r, g, b;
    Color() { r = g = b = 0; }
    Color(lld r, lld g, lld b) { this->r = r, this->g = g, this->b = b; }
    void update(Color &c) { this->r = c.r, this->g = c.g, this->b = c.b; }
    void update(lld r, lld g, lld b) { this->r = r, this->g = g, this->b = b; }
};

class Ray
{
public:
    Point starting_point;
    Vector direction;
    Ray() {}
    Ray(Point &point, Vector &dir) { starting_point.update(point), direction.update(dir); }
    void update(Point &point, Vector &dir) { starting_point.update(point), direction.update(dir); }
};

class BoundingBox
{
public:
    lld length, width, height;
    BoundingBox() { length = width = height = 0; }
    BoundingBox(lld length, lld width, lld height) { this->length = length, this->width = width, this->height = height; }
    void update(BoundingBox &bx) { this->length = bx.length, this->width = bx.width, this->height = bx.height; }
    void update(lld length, lld width, lld height) { this->length = length, this->width = width, this->height = height; }
};

class Coefficients
{
public:
    lld ambient, diffuse, specular, reflection;
    Coefficients() { ambient = diffuse = specular = reflection = 0; }
    Coefficients(lld ambient, lld diffuse, lld specular, lld reflection) { this->ambient = ambient, this->diffuse = diffuse, this->reflection = reflection, this->specular = specular; }
    void update(Coefficients &c) { this->ambient = c.ambient, this->diffuse = c.diffuse, this->reflection = c.reflection, this->specular = c.specular; }
    void update(lld ambient, lld diffuse, lld specular, lld reflection) { this->ambient = ambient, this->diffuse = diffuse, this->reflection = reflection, this->specular = specular; }
};

class PointLight
{
public:
    Point light_pos;
    Color color;

    PointLight() {}
    PointLight(Point &pos, Color &col)
    {
        light_pos.update(pos);
        color.update(col);
    }
};

class SpotLight
{
public:
    PointLight point_light;
    Vector light_direction;
    lld cutoff_angle;

    SpotLight(Point &pos, Color &col, Vector dir, lld angle)
    {
        PointLight pt(pos, col);
        point_light = pt;
        light_direction = dir;
        cutoff_angle = angle;
    }
};

class Object
{
public:
    Point reference_point;     // x, y, z
    Color color;               // RGB
    Coefficients coefficients; // ambient, diffuse, specular, reflection
    int shininess;             // specular exponent

    Object() {};
    virtual void draw() {};
    virtual lld intersect(Ray &ray, Color &color, int level) { return -1.0; }
    void setReferencePoint(Point &po) { reference_point.update(po); }
    void setColor(Color &col) { color.update(col); }
    void setShine(int shine) { shininess = shine; }
    void setCoefficients(Coefficients &co) { coefficients.update(co); }
};

// Global variables
vector<Object *> objects;
vector<PointLight> point_lights;
vector<SpotLight> spot_lights;
int recursion_level;
int image_resolution;

class Sphere : public Object
{
public:
    lld radius;

    Sphere(lld rad, Point &center, Color &color, Coefficients &coefficients, int shininess)
    {
        radius = rad;
        setReferencePoint(center);
        setColor(color);
        setCoefficients(coefficients);
        setShine(shininess);
    }

    void draw() override
    {
        // OpenGL code to draw a sphere centered at reference_Point with given radius
        glPushMatrix();
        glTranslatef(reference_point.x, reference_point.y, reference_point.z);
        glColor3f(color.r, color.g, color.b);
        glutSolidSphere(radius, 50, 50); // Adjust slices and stacks for smoothness
        glPopMatrix();
    }

    // Intersect method for Sphere class
    lld intersect(Ray &ray, Color &col, int level) override
    {
        // Vector from ray origin to sphere center
        Vector oc(ray.starting_point.x - reference_point.x,
                  ray.starting_point.y - reference_point.y,
                  ray.starting_point.z - reference_point.z);

        // Quadratic equation coefficients: at² + bt + c = 0
        lld a = ray.direction.dot(ray.direction);
        lld b = 2.0 * oc.dot(ray.direction);
        lld c = oc.dot(oc) - radius * radius;

        // Calculate discriminant
        lld discriminant = b * b - 4 * a * c;
        if (discriminant < 0)
            return -1.0; // No intersection

        // Calculate both roots
        lld sqrt_discriminant = sqrt(discriminant);
        lld t1 = (-b - sqrt_discriminant) / (2.0 * a);
        lld t2 = (-b + sqrt_discriminant) / (2.0 * a);

        // Choose the nearest positive intersection
        lld t = -1.0;
        if (t1 > EPSILON)
        {
            t = t1;
        }
        else if (t2 > EPSILON)
        {
            t = t2;
        }

        if (t > EPSILON)
        {
            if (level == 0)
            {
                // Set the color output parameter
                col.update(this->color);
                return t;
            }

            // Compute intersection point
            Point intersectionPoint(ray.starting_point.x + t * ray.direction.x,
                                    ray.starting_point.y + t * ray.direction.y,
                                    ray.starting_point.z + t * ray.direction.z);

            // Compute normal
            Vector normal(intersectionPoint.x - reference_point.x,
                          intersectionPoint.y - reference_point.y,
                          intersectionPoint.z - reference_point.z);
            normal = normal.normalize();

            // Ambient component
            col.update(this->color.r * coefficients.ambient,
                       this->color.g * coefficients.ambient,
                       this->color.b * coefficients.ambient);

            // Process point lights
            for (const auto &pl : point_lights)
            {
                Vector lightDir(pl.light_pos.x - intersectionPoint.x,
                                pl.light_pos.y - intersectionPoint.y,
                                pl.light_pos.z - intersectionPoint.z);
                lld distToLight = lightDir.length();
                lightDir = lightDir.normalize();

                Point shadowRayStart(intersectionPoint.x + normal.x * EPSILON,
                                     intersectionPoint.y + normal.y * EPSILON,
                                     intersectionPoint.z + normal.z * EPSILON);
                Ray shadowRay(shadowRayStart, lightDir);

                bool inShadow = false;
                for (auto *obj : objects)
                {
                    Color dummyColor;
                    lld t_shadow = obj->intersect(shadowRay, dummyColor, 0);
                    if (t_shadow > EPSILON && t_shadow < distToLight)
                    {
                        inShadow = true;
                        break;
                    }
                }

                if (!inShadow)
                {
                    // Diffuse
                    lld lambertValue = max((lld)0.0, normal.dot(lightDir));
                    col.r += pl.color.r * coefficients.diffuse * lambertValue * this->color.r;
                    col.g += pl.color.g * coefficients.diffuse * lambertValue * this->color.g;
                    col.b += pl.color.b * coefficients.diffuse * lambertValue * this->color.b;

                    // Specular
                    Vector scaledNormal = normal.scalar_mul(2.0 * normal.dot(lightDir));
                    Vector reflectDir = lightDir.sub(scaledNormal);
                    reflectDir = reflectDir.normalize();
                    Vector viewDir = ray.direction.scalar_mul(-1.0).normalize();
                    lld phongValue = pow(max((lld)0.0, reflectDir.dot(viewDir)), shininess);
                    col.r += pl.color.r * coefficients.specular * phongValue * this->color.r;
                    col.g += pl.color.g * coefficients.specular * phongValue * this->color.g;
                    col.b += pl.color.b * coefficients.specular * phongValue * this->color.b;
                }
            }

            // Process spotlights
            for (const auto &sl : spot_lights)
            {
                Vector lightDir(sl.point_light.light_pos.x - intersectionPoint.x,
                                sl.point_light.light_pos.y - intersectionPoint.y,
                                sl.point_light.light_pos.z - intersectionPoint.z);
                lld distToLight = lightDir.length();
                lightDir = lightDir.normalize();

                Point shadowRayStart(intersectionPoint.x + normal.x * EPSILON,
                                     intersectionPoint.y + normal.y * EPSILON,
                                     intersectionPoint.z + normal.z * EPSILON);
                Ray shadowRay(shadowRayStart, lightDir);

                bool inShadow = false;
                for (auto *obj : objects)
                {
                    Color dummyColor;
                    lld t_shadow = obj->intersect(shadowRay, dummyColor, 0);
                    if (t_shadow > EPSILON && t_shadow < distToLight)
                    {
                        inShadow = true;
                        break;
                    }
                }

                if (!inShadow)
                {
                    // Diffuse
                    lld lambertValue = max((lld)0.0, normal.dot(lightDir));
                    col.r += sl.point_light.color.r * coefficients.diffuse * lambertValue * this->color.r;
                    col.g += sl.point_light.color.g * coefficients.diffuse * lambertValue * this->color.g;
                    col.b += sl.point_light.color.b * coefficients.diffuse * lambertValue * this->color.b;

                    // Specular
                    Vector scaledNormal = normal.scalar_mul(2.0 * normal.dot(lightDir));
                    Vector reflectDir = lightDir.sub(scaledNormal);
                    reflectDir = reflectDir.normalize();
                    Vector viewDir = ray.direction.scalar_mul(-1.0).normalize();
                    lld phongValue = pow(max((lld)0.0, reflectDir.dot(viewDir)), shininess);
                    col.r += sl.point_light.color.r * coefficients.specular * phongValue * this->color.r;
                    col.g += sl.point_light.color.g * coefficients.specular * phongValue * this->color.g;
                    col.b += sl.point_light.color.b * coefficients.specular * phongValue * this->color.b;
                }
            }

            // Recursive reflection - optimized version
            if (level < recursion_level && coefficients.reflection > EPSILON)
            {
                Vector viewDir = ray.direction.scalar_mul(-1.0).normalize();
                Vector scaledNormal = normal.scalar_mul(2.0 * normal.dot(ray.direction));
                Vector reflectDir = ray.direction.sub(scaledNormal).normalize();
                Point reflectStart(intersectionPoint.x + reflectDir.x * EPSILON,
                                   intersectionPoint.y + reflectDir.y * EPSILON,
                                   intersectionPoint.z + reflectDir.z * EPSILON);
                Ray reflectRay(reflectStart, reflectDir);
                Color reflectColor(0, 0, 0);

                // Find nearest intersecting object
                Object *nearestObj = nullptr;
                lld tMin = 1e9;
                for (auto *obj : objects)
                {
                    Color objColor;
                    lld t_reflect = obj->intersect(reflectRay, objColor, 0); // Only check intersection, no lighting
                    if (t_reflect > EPSILON && t_reflect < tMin)
                    {
                        tMin = t_reflect;
                        nearestObj = obj;
                    }
                }

                // Only compute lighting if we hit something and reflection is significant
                if (nearestObj && coefficients.reflection > 0.1)
                {
                    Color finalReflectColor;
                    nearestObj->intersect(reflectRay, finalReflectColor, level + 1);

                    // Scale by reflection coefficient
                    col.r += finalReflectColor.r * coefficients.reflection;
                    col.g += finalReflectColor.g * coefficients.reflection;
                    col.b += finalReflectColor.b * coefficients.reflection;
                }
            }

            // Clamp colors
            col.r = min((lld)1.0, max((lld)0.0, col.r));
            col.g = min((lld)1.0, max((lld)0.0, col.g));
            col.b = min((lld)1.0, max((lld)0.0, col.b));

            return t;
        }

        return -1.0; // No valid intersection
    }
};

class Triangle : public Object
{
public:
    Point p1, p2, p3;

    Triangle(Point &p1, Point &p2, Point &p3, Color &color, Coefficients &coefficients, int shininess)
    {
        this->p1.update(p1), this->p2.update(p2), this->p3.update(p3);
        setReferencePoint(p1);
        setColor(color);
        setCoefficients(coefficients);
        setShine(shininess);
    }

    void draw() override
    {
        glPushMatrix();
        glColor3f(color.r, color.g, color.b);
        glBegin(GL_TRIANGLES);
        glVertex3f(p1.x, p1.y, p1.z);
        glVertex3f(p2.x, p2.y, p2.z);
        glVertex3f(p3.x, p3.y, p3.z);
        glEnd();
        glPopMatrix();
    }

    // Triangle intersect function with Phong lighting and recursive reflection
    lld intersect(Ray &ray, Color &col, int level) override
    {
        // Edge vectors
        Vector edge1(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
        Vector edge2(p3.x - p1.x, p3.y - p1.y, p3.z - p1.z);

        // Begin calculating determinant - also used to calculate u parameter
        Vector h = ray.direction.cross(edge2);
        lld det = edge1.dot(h);

        // If determinant is near zero, ray lies in plane of triangle
        if (det > -EPSILON && det < EPSILON)
        {
            return -1.0;
        }

        lld inv_det = 1.0 / det;

        // Calculate distance from vertex p1 to ray origin
        Vector s(ray.starting_point.x - p1.x, ray.starting_point.y - p1.y, ray.starting_point.z - p1.z);

        // Calculate u parameter and test bounds
        lld u = s.dot(h) * inv_det;
        if (u < 0.0 || u > 1.0)
        {
            return -1.0;
        }

        // Prepare to test v parameter
        Vector q = s.cross(edge1);

        // Calculate v parameter and test bounds
        lld v = ray.direction.dot(q) * inv_det;
        if (v < 0.0 || u + v > 1.0)
        {
            return -1.0;
        }

        // Calculate t parameter (distance along ray)
        lld t = edge2.dot(q) * inv_det;

        if (t > EPSILON)
        {
            if (level == 0)
            {
                col.update(this->color);
                return t;
            }

            // Phong lighting model implementation
            Point intersectionPoint(ray.starting_point.x + t * ray.direction.x,
                                    ray.starting_point.y + t * ray.direction.y,
                                    ray.starting_point.z + t * ray.direction.z);

            // Calculate normal using cross product of edge vectors
            Vector normal = edge1.cross(edge2).normalize();

            // Ambient component
            col.update(this->color.r * coefficients.ambient,
                       this->color.g * coefficients.ambient,
                       this->color.b * coefficients.ambient);

            // Process point lights
            for (const auto &pl : point_lights)
            {
                Vector lightDir(pl.light_pos.x - intersectionPoint.x,
                                pl.light_pos.y - intersectionPoint.y,
                                pl.light_pos.z - intersectionPoint.z);
                lld distToLight = lightDir.length();
                lightDir = lightDir.normalize();

                Point shadowRayStart(intersectionPoint.x + normal.x * EPSILON,
                                     intersectionPoint.y + normal.y * EPSILON,
                                     intersectionPoint.z + normal.z * EPSILON);
                Ray shadowRay(shadowRayStart, lightDir);

                bool inShadow = false;
                for (auto *obj : objects)
                {
                    Color dummyColor;
                    lld t_shadow = obj->intersect(shadowRay, dummyColor, 0);
                    if (t_shadow > EPSILON && t_shadow < distToLight)
                    {
                        inShadow = true;
                        break;
                    }
                }

                if (!inShadow)
                {
                    // Diffuse
                    lld lambertValue = max((lld)0.0, normal.dot(lightDir));
                    col.r += pl.color.r * coefficients.diffuse * lambertValue * this->color.r;
                    col.g += pl.color.g * coefficients.diffuse * lambertValue * this->color.g;
                    col.b += pl.color.b * coefficients.diffuse * lambertValue * this->color.b;

                    // Specular
                    Vector scaledNormal = normal.scalar_mul(2.0 * normal.dot(lightDir));
                    Vector reflectDir = lightDir.sub(scaledNormal);
                    reflectDir = reflectDir.normalize();
                    Vector viewDir = ray.direction.scalar_mul(-1.0).normalize();
                    lld phongValue = pow(max((lld)0.0, reflectDir.dot(viewDir)), shininess);
                    col.r += pl.color.r * coefficients.specular * phongValue * this->color.r;
                    col.g += pl.color.g * coefficients.specular * phongValue * this->color.g;
                    col.b += pl.color.b * coefficients.specular * phongValue * this->color.b;
                }
            }

            // Process spotlights
            for (const auto &sl : spot_lights)
            {
                Vector lightDir(sl.point_light.light_pos.x - intersectionPoint.x,
                                sl.point_light.light_pos.y - intersectionPoint.y,
                                sl.point_light.light_pos.z - intersectionPoint.z);
                lld distToLight = lightDir.length();
                lightDir = lightDir.normalize();

                Point shadowRayStart(intersectionPoint.x + normal.x * EPSILON,
                                     intersectionPoint.y + normal.y * EPSILON,
                                     intersectionPoint.z + normal.z * EPSILON);
                Ray shadowRay(shadowRayStart, lightDir);

                bool inShadow = false;
                for (auto *obj : objects)
                {
                    Color dummyColor;
                    lld t_shadow = obj->intersect(shadowRay, dummyColor, 0);
                    if (t_shadow > EPSILON && t_shadow < distToLight)
                    {
                        inShadow = true;
                        break;
                    }
                }

                if (!inShadow)
                {
                    // Diffuse
                    lld lambertValue = max((lld)0.0, normal.dot(lightDir));
                    col.r += sl.point_light.color.r * coefficients.diffuse * lambertValue * this->color.r;
                    col.g += sl.point_light.color.g * coefficients.diffuse * lambertValue * this->color.g;
                    col.b += sl.point_light.color.b * coefficients.diffuse * lambertValue * this->color.b;

                    // Specular
                    Vector scaledNormal = normal.scalar_mul(2.0 * normal.dot(lightDir));
                    Vector reflectDir = lightDir.sub(scaledNormal);
                    reflectDir = reflectDir.normalize();
                    Vector viewDir = ray.direction.scalar_mul(-1.0).normalize();
                    lld phongValue = pow(max((lld)0.0, reflectDir.dot(viewDir)), shininess);
                    col.r += sl.point_light.color.r * coefficients.specular * phongValue * this->color.r;
                    col.g += sl.point_light.color.g * coefficients.specular * phongValue * this->color.g;
                    col.b += sl.point_light.color.b * coefficients.specular * phongValue * this->color.b;
                }
            }

            // Recursive reflection - optimized version
            if (level < recursion_level && coefficients.reflection > EPSILON)
            {
                Vector viewDir = ray.direction.scalar_mul(-1.0).normalize();
                Vector scaledNormal = normal.scalar_mul(2.0 * normal.dot(ray.direction));
                Vector reflectDir = ray.direction.sub(scaledNormal).normalize();
                Point reflectStart(intersectionPoint.x + reflectDir.x * EPSILON,
                                   intersectionPoint.y + reflectDir.y * EPSILON,
                                   intersectionPoint.z + reflectDir.z * EPSILON);
                Ray reflectRay(reflectStart, reflectDir);
                Color reflectColor(0, 0, 0);

                // Find nearest intersecting object
                Object *nearestObj = nullptr;
                lld tMin = 1e9;
                for (auto *obj : objects)
                {
                    Color objColor;
                    lld t_reflect = obj->intersect(reflectRay, objColor, 0); // Only check intersection, no lighting
                    if (t_reflect > EPSILON && t_reflect < tMin)
                    {
                        tMin = t_reflect;
                        nearestObj = obj;
                    }
                }

                // Only compute lighting if we hit something and reflection is significant
                if (nearestObj && coefficients.reflection > 0.1)
                {
                    Color finalReflectColor;
                    nearestObj->intersect(reflectRay, finalReflectColor, level + 1);

                    // Scale by reflection coefficient
                    col.r += finalReflectColor.r * coefficients.reflection;
                    col.g += finalReflectColor.g * coefficients.reflection;
                    col.b += finalReflectColor.b * coefficients.reflection;
                }
            }

            // Clamp colors
            col.r = min((lld)1.0, max((lld)0.0, col.r));
            col.g = min((lld)1.0, max((lld)0.0, col.g));
            col.b = min((lld)1.0, max((lld)0.0, col.b));

            return t;
        }

        return -1.0; // Ray intersection, but behind ray origin
    }
};

class GeneralQuadric : public Object
{
public:
    lld A, B, C, D, E, F, G, H, I, J; // Coefficients of quadric equation
    BoundingBox bounding_box;         // Bounding box dimensions
    Point bounding_box_reference;     // Reference point for bounding box

    GeneralQuadric(lld A, lld B, lld C, lld D, lld E, lld F, lld G, lld H, lld I, lld J,
                   Point &box_reference, lld box_length, lld box_width, lld box_height, Color &color, Coefficients &coefficients, int shininess)
    {
        this->A = A, this->B = B, this->C = C, this->D = D, this->E = E, this->F = F, this->G = G, this->H = H, this->I = I, this->J = J;
        bounding_box_reference.update(box_reference);
        bounding_box.update(box_length, box_width, box_height);
        setReferencePoint(box_reference);
        setColor(color);
        setCoefficients(coefficients);
        setShine(shininess);
    }

    void draw() override
    {
        // General quadric surfaces are not drawn in OpenGL; only for ray tracing
    }

    // GeneralQuadric intersect function with Phong lighting and recursive reflection
    lld intersect(Ray &ray, Color &col, int level) override
    {
        // Ray origin and direction
        lld ox = ray.starting_point.x;
        lld oy = ray.starting_point.y;
        lld oz = ray.starting_point.z;
        lld dx = ray.direction.x;
        lld dy = ray.direction.y;
        lld dz = ray.direction.z;

        // Solve the quadric equation
        // General quadric: Ax² + By² + Cz² + Dxy + Exz + Fyz + Gx + Hy + Iz + J = 0
        // Substitute ray equation: P = O + t*D where O is origin, D is direction

        // Quadratic equation coefficients: at² + bt + c = 0
        lld a = A * dx * dx + B * dy * dy + C * dz * dz +
                D * dx * dy + E * dx * dz + F * dy * dz;

        lld b = 2.0 * A * ox * dx + 2.0 * B * oy * dy + 2.0 * C * oz * dz +
                D * (ox * dy + oy * dx) + E * (ox * dz + oz * dx) + F * (oy * dz + oz * dy) +
                G * dx + H * dy + I * dz;

        lld c = A * ox * ox + B * oy * oy + C * oz * oz +
                D * ox * oy + E * ox * oz + F * oy * oz +
                G * ox + H * oy + I * oz + J;

        // Handle different cases based on coefficient 'a'
        vector<lld> valid_t_values;

        if (abs(a) < EPSILON)
        {
            // Linear equation: bt + c = 0
            if (abs(b) < EPSILON)
            {
                return -1.0; // No solution or infinitely many solutions
            }
            lld t = -c / b;
            if (t > EPSILON)
            {
                valid_t_values.push_back(t);
            }
        }
        else
        {
            // Quadratic equation
            lld discriminant = b * b - 4.0 * a * c;

            if (discriminant < 0)
            {
                return -1.0; // No real solutions
            }

            lld sqrt_discriminant = sqrt(discriminant);
            lld t1 = (-b - sqrt_discriminant) / (2.0 * a);
            lld t2 = (-b + sqrt_discriminant) / (2.0 * a);

            if (t1 > EPSILON)
                valid_t_values.push_back(t1);
            if (t2 > EPSILON)
                valid_t_values.push_back(t2);
        }

        // Check each valid t value against bounding box constraints
        for (lld t : valid_t_values)
        {
            // Calculate intersection point
            lld px = ox + t * dx;
            lld py = oy + t * dy;
            lld pz = oz + t * dz;

            // Check bounding box constraints
            bool within_bounds = true;

            // Check x-bounds only if length > 0 (non-zero means clipping is active)
            if (bounding_box.length > EPSILON)
            {
                lld xmin = bounding_box_reference.x;
                lld xmax = bounding_box_reference.x + bounding_box.length;
                if (px < xmin || px > xmax)
                {
                    within_bounds = false;
                }
            }

            // Check y-bounds only if width > 0 (non-zero means clipping is active)
            if (within_bounds && bounding_box.width > EPSILON)
            {
                lld ymin = bounding_box_reference.y;
                lld ymax = bounding_box_reference.y + bounding_box.width;
                if (py < ymin || py > ymax)
                {
                    within_bounds = false;
                }
            }

            // Check z-bounds only if height > 0 (non-zero means clipping is active)
            if (within_bounds && bounding_box.height > EPSILON)
            {
                lld zmin = bounding_box_reference.z;
                lld zmax = bounding_box_reference.z + bounding_box.height;
                if (pz < zmin || pz > zmax)
                {
                    within_bounds = false;
                }
            }

            if (within_bounds)
            {
                if (level == 0)
                {
                    col.update(this->color);
                    return t;
                }

                // Phong lighting model implementation
                Point intersectionPoint(px, py, pz);

                // Calculate normal using gradient of quadric surface
                // ∇F = (2Ax + Dy + Ez + G, 2By + Dx + Fz + H, 2Cz + Ex + Fy + I)
                Vector normal(2.0 * A * px + D * py + E * pz + G,
                              2.0 * B * py + D * px + F * pz + H,
                              2.0 * C * pz + E * px + F * py + I);
                normal = normal.normalize();

                // Ambient component
                col.update(this->color.r * coefficients.ambient,
                           this->color.g * coefficients.ambient,
                           this->color.b * coefficients.ambient);

                // Process point lights
                for (const auto &pl : point_lights)
                {
                    Vector lightDir(pl.light_pos.x - intersectionPoint.x,
                                    pl.light_pos.y - intersectionPoint.y,
                                    pl.light_pos.z - intersectionPoint.z);
                    lld distToLight = lightDir.length();
                    lightDir = lightDir.normalize();

                    Point shadowRayStart(intersectionPoint.x + normal.x * EPSILON,
                                         intersectionPoint.y + normal.y * EPSILON,
                                         intersectionPoint.z + normal.z * EPSILON);
                    Ray shadowRay(shadowRayStart, lightDir);

                    bool inShadow = false;
                    for (auto *obj : objects)
                    {
                        Color dummyColor;
                        lld t_shadow = obj->intersect(shadowRay, dummyColor, 0);
                        if (t_shadow > EPSILON && t_shadow < distToLight)
                        {
                            inShadow = true;
                            break;
                        }
                    }

                    if (!inShadow)
                    {
                        // Diffuse
                        lld lambertValue = max((lld)0.0, normal.dot(lightDir));
                        col.r += pl.color.r * coefficients.diffuse * lambertValue * this->color.r;
                        col.g += pl.color.g * coefficients.diffuse * lambertValue * this->color.g;
                        col.b += pl.color.b * coefficients.diffuse * lambertValue * this->color.b;

                        // Specular
                        Vector scaledNormal = normal.scalar_mul(2.0 * normal.dot(lightDir));
                        Vector reflectDir = lightDir.sub(scaledNormal);
                        reflectDir = reflectDir.normalize();
                        Vector viewDir = ray.direction.scalar_mul(-1.0).normalize();
                        lld phongValue = pow(max((lld)0.0, reflectDir.dot(viewDir)), shininess);
                        col.r += pl.color.r * coefficients.specular * phongValue * this->color.r;
                        col.g += pl.color.g * coefficients.specular * phongValue * this->color.g;
                        col.b += pl.color.b * coefficients.specular * phongValue * this->color.b;
                    }
                }

                // Process spotlights
                for (const auto &sl : spot_lights)
                {
                    Vector lightDir(sl.point_light.light_pos.x - intersectionPoint.x,
                                    sl.point_light.light_pos.y - intersectionPoint.y,
                                    sl.point_light.light_pos.z - intersectionPoint.z);
                    lld distToLight = lightDir.length();
                    lightDir = lightDir.normalize();

                    Point shadowRayStart(intersectionPoint.x + normal.x * EPSILON,
                                         intersectionPoint.y + normal.y * EPSILON,
                                         intersectionPoint.z + normal.z * EPSILON);
                    Ray shadowRay(shadowRayStart, lightDir);

                    bool inShadow = false;
                    for (auto *obj : objects)
                    {
                        Color dummyColor;
                        lld t_shadow = obj->intersect(shadowRay, dummyColor, 0);
                        if (t_shadow > EPSILON && t_shadow < distToLight)
                        {
                            inShadow = true;
                            break;
                        }
                    }

                    if (!inShadow)
                    {
                        // Diffuse
                        lld lambertValue = max((lld)0.0, normal.dot(lightDir));
                        col.r += sl.point_light.color.r * coefficients.diffuse * lambertValue * this->color.r;
                        col.g += sl.point_light.color.g * coefficients.diffuse * lambertValue * this->color.g;
                        col.b += sl.point_light.color.b * coefficients.diffuse * lambertValue * this->color.b;

                        // Specular
                        Vector scaledNormal = normal.scalar_mul(2.0 * normal.dot(lightDir));
                        Vector reflectDir = lightDir.sub(scaledNormal);
                        reflectDir = reflectDir.normalize();
                        Vector viewDir = ray.direction.scalar_mul(-1.0).normalize();
                        lld phongValue = pow(max((lld)0.0, reflectDir.dot(viewDir)), shininess);
                        col.r += sl.point_light.color.r * coefficients.specular * phongValue * this->color.r;
                        col.g += sl.point_light.color.g * coefficients.specular * phongValue * this->color.g;
                        col.b += sl.point_light.color.b * coefficients.specular * phongValue * this->color.b;
                    }
                }

                // Recursive reflection - optimized version
                if (level < recursion_level && coefficients.reflection > EPSILON)
                {
                    Vector viewDir = ray.direction.scalar_mul(-1.0).normalize();
                    Vector scaledNormal = normal.scalar_mul(2.0 * normal.dot(ray.direction));
                    Vector reflectDir = ray.direction.sub(scaledNormal).normalize();
                    Point reflectStart(intersectionPoint.x + reflectDir.x * EPSILON,
                                       intersectionPoint.y + reflectDir.y * EPSILON,
                                       intersectionPoint.z + reflectDir.z * EPSILON);
                    Ray reflectRay(reflectStart, reflectDir);
                    Color reflectColor(0, 0, 0);

                    // Find nearest intersecting object
                    Object *nearestObj = nullptr;
                    lld tMin = 1e9;
                    for (auto *obj : objects)
                    {
                        Color objColor;
                        lld t_reflect = obj->intersect(reflectRay, objColor, 0); // Only check intersection, no lighting
                        if (t_reflect > EPSILON && t_reflect < tMin)
                        {
                            tMin = t_reflect;
                            nearestObj = obj;
                        }
                    }

                    // Only compute lighting if we hit something and reflection is significant
                    if (nearestObj && coefficients.reflection > 0.1)
                    {
                        Color finalReflectColor;
                        nearestObj->intersect(reflectRay, finalReflectColor, level + 1);

                        // Scale by reflection coefficient
                        col.r += finalReflectColor.r * coefficients.reflection;
                        col.g += finalReflectColor.g * coefficients.reflection;
                        col.b += finalReflectColor.b * coefficients.reflection;
                    }
                }

                // Clamp colors
                col.r = min((lld)1.0, max((lld)0.0, col.r));
                col.g = min((lld)1.0, max((lld)0.0, col.g));
                col.b = min((lld)1.0, max((lld)0.0, col.b));

                return t;
            }
        }

        return -1.0; // No valid intersection within bounding constraints
    }
};

class Floor : public Object
{
public:
    lld floorWidth, tileWidth;
    Color color1, color2; // Alternating colors for checkerboard

    // Texture-related members
    bool useTexture;
    unsigned char *textureData;
    int textureWidth, textureHeight, textureChannels;
    string texturePath;

    Floor(lld fWidth, lld tWidth, Coefficients &co)
    {
        floorWidth = fWidth;
        tileWidth = tWidth;
        reference_point.update(-fWidth / 2, -fWidth / 2, 0);
        color1.update(1.0, 1.0, 1.0); // White
        color2.update(0.0, 0.0, 0.0); // Black
        setCoefficients(co);          // Example coefficients
        setShine(5);                  // Example shine

        // Initialize texture-related members
        useTexture = true;
        textureData = nullptr;
        textureWidth = textureHeight = textureChannels = 0;
    }

    bool loadTexture(const string &path)
    {
        texturePath = path;
        textureData = stbi_load(path.c_str(), &textureWidth, &textureHeight, &textureChannels, 0);

        if (textureData == nullptr)
        {
            cout << "Failed to load texture: " << path << endl;
            useTexture = false;
            return false;
        }

        useTexture = true;
        cout << "Loaded texture" << endl;
        return true;
    }

    Color getTextureColor(lld x, lld y)
    {
        if (!useTexture || textureData == nullptr)
        {
            // Return checkerboard pattern if no texture
            int tileX = (int)((x - reference_point.x) / tileWidth);
            int tileY = (int)((y - reference_point.y) / tileWidth);
            return ((tileX + tileY) % 2 == 0) ? color1 : color2;
        }

        // Convert world coordinates to texture coordinates
        lld u = (x - reference_point.x) / floorWidth;
        lld v = (y - reference_point.y) / floorWidth;

        u = max(0.0, min(1.0, u));
        v = max(0.0, min(1.0, v));

        int px = (int)(u * (textureWidth - 1));
        int py = (int)(v * (textureHeight - 1));

        px = max(0, min(textureWidth - 1, px));
        py = max(0, min(textureHeight - 1, py));

        int index = (py * textureWidth + px) * textureChannels;

        if (textureChannels >= 3)
        {
            lld r = textureData[index] / 255.0;
            lld g = textureData[index + 1] / 255.0;
            lld b = textureData[index + 2] / 255.0;
            return Color(r, g, b);
        }
        else if (textureChannels == 1)
        {
            lld gray = textureData[index] / 255.0;
            return Color(gray, gray, gray);
        }

        return Color(1.0, 1.0, 1.0);
    }

    ~Floor()
    {
        if (textureData != nullptr)
        {
            stbi_image_free(textureData);
        }
    }

    void draw() override
    {
        glPushMatrix();
        glTranslatef(reference_point.x, reference_point.y, reference_point.z);
        for (lld x = 0; x < floorWidth; x += tileWidth)
        {
            for (lld y = 0; y < floorWidth; y += tileWidth)
            {
                // Determine color based on tile position
                if (((int)(x / tileWidth) + (int)(y / tileWidth)) % 2 == 0)
                {
                    glColor3f(color1.r, color1.g, color1.b);
                }
                else
                {
                    glColor3f(color2.r, color2.g, color2.b);
                }
                glBegin(GL_QUADS);
                glVertex3f(x, y, 0);
                glVertex3f(x + tileWidth, y, 0);
                glVertex3f(x + tileWidth, y + tileWidth, 0);
                glVertex3f(x, y + tileWidth, 0);
                glEnd();
            }
        }
        glPopMatrix();
    }

    // Floor intersect function with Phong lighting and recursive reflection
    lld intersect(Ray &ray, Color &col, int level) override
    {
        // Floor is at z = 0, normal is (0, 0, 1)
        // Ray equation: P = origin + t * direction
        // Plane equation: z = 0

        if (abs(ray.direction.z) < EPSILON)
        {
            return -1.0; // Ray is parallel to floor
        }

        // t = (0 - ray.starting_point.z) / ray.direction.z
        lld t = -ray.starting_point.z / ray.direction.z;

        if (t < EPSILON)
        {
            return -1.0; // Intersection behind ray origin
        }

        // Calculate intersection point
        lld x = ray.starting_point.x + t * ray.direction.x;
        lld y = ray.starting_point.y + t * ray.direction.y;

        // Check if intersection is within floor bounds
        if (x < reference_point.x || x > reference_point.x + floorWidth ||
            y < reference_point.y || y > reference_point.y + floorWidth)
        {
            return -1.0; // Outside floor bounds
        }

        Color tileColor = getTextureColor(x, y);

        if (level == 0)
        {
            col.update(tileColor);
            return t;
        }

        // Phong lighting model implementation
        Point intersectionPoint(x, y, 0);
        Vector normal(0, 0, 1); // Floor normal is always (0, 0, 1)

        // Ambient component
        col.update(tileColor.r * coefficients.ambient,
                   tileColor.g * coefficients.ambient,
                   tileColor.b * coefficients.ambient);

        // Process point lights
        for (const auto &pl : point_lights)
        {
            Vector lightDir(pl.light_pos.x - intersectionPoint.x,
                            pl.light_pos.y - intersectionPoint.y,
                            pl.light_pos.z - intersectionPoint.z);
            lld distToLight = lightDir.length();
            lightDir = lightDir.normalize();

            Ray shadowRay(intersectionPoint, lightDir);

            bool inShadow = false;
            for (auto *obj : objects)
            {
                Color dummyColor;
                lld t_shadow = obj->intersect(shadowRay, dummyColor, 0);
                if (t_shadow > EPSILON && t_shadow < distToLight)
                {
                    inShadow = true;
                    break;
                }
            }

            if (!inShadow)
            {
                // Diffuse
                lld lambertValue = max((lld)0.0, normal.dot(lightDir));
                col.r += pl.color.r * coefficients.diffuse * lambertValue * tileColor.r;
                col.g += pl.color.g * coefficients.diffuse * lambertValue * tileColor.g;
                col.b += pl.color.b * coefficients.diffuse * lambertValue * tileColor.b;

                // Specular
                Vector scaledNormal = normal.scalar_mul(2.0 * normal.dot(lightDir));
                Vector reflectDir = lightDir.sub(scaledNormal);
                reflectDir = reflectDir.normalize();
                Vector viewDir = ray.direction.scalar_mul(-1.0).normalize();
                lld phongValue = pow(max((lld)0.0, reflectDir.dot(viewDir)), shininess);
                col.r += pl.color.r * coefficients.specular * phongValue * tileColor.r;
                col.g += pl.color.g * coefficients.specular * phongValue * tileColor.g;
                col.b += pl.color.b * coefficients.specular * phongValue * tileColor.b;
            }
        }

        // Process spotlights
        for (const auto &sl : spot_lights)
        {
            Vector lightDir(sl.point_light.light_pos.x - intersectionPoint.x,
                            sl.point_light.light_pos.y - intersectionPoint.y,
                            sl.point_light.light_pos.z - intersectionPoint.z);
            lld distToLight = lightDir.length();
            lightDir = lightDir.normalize();

            Vector spotDir = sl.light_direction;
            spotDir.update(spotDir.normalize());
            Vector negSpotDir = spotDir.scalar_mul(-1.0);
            lld cosAngle = lightDir.dot(negSpotDir);
            if (cosAngle < cos(sl.cutoff_angle * M_PI / 180.0))
                continue;

            Ray shadowRay(intersectionPoint, lightDir);
            bool inShadow = false;
            for (auto *obj : objects)
            {
                Color dummyColor;
                lld t_shadow = obj->intersect(shadowRay, dummyColor, 0);
                if (t_shadow > EPSILON && t_shadow < distToLight)
                {
                    inShadow = true;
                    break;
                }
            }

            if (!inShadow)
            {
                // Diffuse
                lld lambertValue = max((lld)0.0, normal.dot(lightDir));
                col.r += sl.point_light.color.r * coefficients.diffuse * lambertValue * tileColor.r;
                col.g += sl.point_light.color.g * coefficients.diffuse * lambertValue * tileColor.g;
                col.b += sl.point_light.color.b * coefficients.diffuse * lambertValue * tileColor.b;

                // Specular
                Vector scaledNormal = normal.scalar_mul(2.0 * normal.dot(lightDir));
                Vector reflectDir = lightDir.sub(scaledNormal);
                reflectDir = reflectDir.normalize();
                Vector viewDir = ray.direction.scalar_mul(-1.0).normalize();
                lld phongValue = pow(max((lld)0.0, reflectDir.dot(viewDir)), shininess);
                col.r += sl.point_light.color.r * coefficients.specular * phongValue * tileColor.r;
                col.g += sl.point_light.color.g * coefficients.specular * phongValue * tileColor.g;
                col.b += sl.point_light.color.b * coefficients.specular * phongValue * tileColor.b;
            }
        }

        // Recursive reflection - optimized version
        if (level < recursion_level && coefficients.reflection > EPSILON)
        {
            Vector viewDir = ray.direction.scalar_mul(-1.0).normalize();
            Vector scaledNormal = normal.scalar_mul(2.0 * normal.dot(ray.direction));
            Vector reflectDir = ray.direction.sub(scaledNormal).normalize();
            Point reflectStart(intersectionPoint.x + reflectDir.x * EPSILON,
                               intersectionPoint.y + reflectDir.y * EPSILON,
                               intersectionPoint.z + reflectDir.z * EPSILON);
            Ray reflectRay(reflectStart, reflectDir);
            Color reflectColor(0, 0, 0);

            // Find nearest intersecting object
            Object *nearestObj = nullptr;
            lld tMin = 1e9;
            for (auto *obj : objects)
            {
                Color objColor;
                lld t_reflect = obj->intersect(reflectRay, objColor, 0); // Only check intersection, no lighting
                if (t_reflect > EPSILON && t_reflect < tMin)
                {
                    tMin = t_reflect;
                    nearestObj = obj;
                }
            }

            // Only compute lighting if we hit something and reflection is significant
            if (nearestObj && coefficients.reflection > 0.1)
            {
                Color finalReflectColor;
                nearestObj->intersect(reflectRay, finalReflectColor, level + 1);

                // Scale by reflection coefficient
                col.r += finalReflectColor.r * coefficients.reflection;
                col.g += finalReflectColor.g * coefficients.reflection;
                col.b += finalReflectColor.b * coefficients.reflection;
            }
        }

        // Clamp colors
        col.r = min((lld)1.0, max((lld)0.0, col.r));
        col.g = min((lld)1.0, max((lld)0.0, col.g));
        col.b = min((lld)1.0, max((lld)0.0, col.b));

        return t;
    }
};