#include "bitmap_image.hpp"
#include <bits/stdc++.h>
using namespace std;
#define endl '\n'
#define gap ' '
typedef long long ll;
typedef long double lld;
const ll infinite = INT64_MAX;

lld eyeX, eyeY, eyeZ;
lld lookX, lookY, lookZ;
lld upX, upY, upZ;
lld fovY, aspectRatio, near, far;
int Screen_Width, Screen_Height;
lld left_limit_of_X, right_limit_of_X;
lld bottom_limit_of_Y, top_limit_of_Y;
lld z_min, z_max;

// for random color
int _random()
{
    return (rand()) % 256;
}

// areat of a triangle given 3 points (2D)
lld area_of_triangle(lld x1, lld y1, lld x2, lld y2, lld x3, lld y3)
{
    return 0.5 * abs(x1 * y2 + x2 * y3 + x3 * y1 - y1 * x2 - y2 * x3 - y3 * x1);
}

// cross product
vector<lld> cross(vector<lld> a, vector<lld> b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// matrix multiplication
vector<vector<lld>> multiply(vector<vector<lld>> A, vector<vector<lld>> B)
{
    vector<vector<lld>> result(4, vector<lld>(4, 0));
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

// discarding small values
lld discard(lld value)
{
    lld integerPart = floorf64x(value);
    lld fractionPart = value - integerPart;
    if (abs(fractionPart) < 1e-8)
    {
        return integerPart;
    }
    return value;
}

// given a point and a transformation matrix it returns the transformed point
vector<lld> transform_a_point(lld px, lld py, lld pz, vector<vector<lld>> transformation_matrix)
{
    vector<lld> point = {px, py, pz, 1};
    vector<lld> result = {0, 0, 0, 0};
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            result[i] += transformation_matrix[i][j] * point[j];
        }
    }
    return {discard(result[0] / result[3]), discard(result[1] / result[3]), discard(result[2] / result[3])};
}

// given a triangle and a point it returns the z value
lld z_value(vector<vector<lld>> triangle, pair<lld, lld> point)
{
    lld area_of_PAB = area_of_triangle(point.first, point.second, triangle[0][0], triangle[0][1], triangle[1][0], triangle[1][1]);
    lld area_of_PBC = area_of_triangle(point.first, point.second, triangle[1][0], triangle[1][1], triangle[2][0], triangle[2][1]);
    lld area_of_PCA = area_of_triangle(point.first, point.second, triangle[2][0], triangle[2][1], triangle[0][0], triangle[0][1]);
    lld area_of_ABC = area_of_triangle(triangle[0][0], triangle[0][1], triangle[1][0], triangle[1][1], triangle[2][0], triangle[2][1]);
    lld diff = area_of_ABC - (area_of_PAB + area_of_PBC + area_of_PCA), z_val = infinite;
    // diff is for checking if the point is inside the triangel in the 2D XY plane
    if (abs(diff) < 1e-8)
    {
        // z value calculated using Barycentric Coordinates system
        lld alpha = area_of_PBC / area_of_ABC;
        lld beta = area_of_PCA / area_of_ABC;
        lld gamma = area_of_PAB / area_of_ABC;
        z_val = alpha * triangle[0][2] + beta * triangle[1][2] + gamma * triangle[2][2];
    }
    return z_val;
}

// data structure for mataining stack
class transformation_stack
{
private:
    vector<vector<vector<lld>>> v_stack{1, vector<vector<lld>>(4, vector<lld>(4, 0))};
    vector<vector<lld>> identity_matrix{4, vector<lld>(4, 0)};

    // Rodrigues Formulla
    vector<lld> R(vector<lld> x, vector<lld> a, lld angle)
    {
        lld a_dot_x = a[0] * x[0] + a[1] * x[1] + a[2] * x[2];
        vector<lld> a_cross_x = cross(a, x);
        vector<lld> rot_v = {0, 0, 0};
        lld theta = angle * M_PI / 180L;
        for (int i = 0; i < 3; i++)
        {
            rot_v[i] = cosl(theta) * x[i] + (1 - cosl(theta)) * a_dot_x * a[i] + sinl(theta) * a_cross_x[i];
        }
        return rot_v;
    }

public:
    transformation_stack()
    {
        // initialized with identity matrix
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                if (i == j)
                {
                    identity_matrix[i][j] = 1;
                }
            }
        }
        v_stack[0] = identity_matrix;
    }
    void push()
    {
        v_stack.push_back(v_stack.back());
    }
    void pop()
    {
        v_stack.pop_back();
    }
    // calculating transalation matrix
    void translate(lld tx, lld ty, lld tz)
    {
        vector<vector<lld>> translation_matrix = identity_matrix;
        translation_matrix[0][3] = tx;
        translation_matrix[1][3] = ty;
        translation_matrix[2][3] = tz;
        v_stack[v_stack.size() - 1] = multiply(v_stack.back(), translation_matrix);
    }
    // calculating scaling matrix
    void scale(lld sx, lld sy, lld sz)
    {
        vector<vector<lld>> scaling_matrix = identity_matrix;
        scaling_matrix[0][0] = sx;
        scaling_matrix[1][1] = sy;
        scaling_matrix[2][2] = sz;
        v_stack[v_stack.size() - 1] = multiply(v_stack.back(), scaling_matrix);
    }
    // calculating rotation matrix
    void rotate(lld angle, lld ax, lld ay, lld az)
    {
        lld l = sqrt(ax * ax + ay * ay + az * az);
        ax /= l, ay /= l, az /= l;
        vector<lld> ci = R({1, 0, 0}, {ax, ay, az}, angle);
        vector<lld> cj = R({0, 1, 0}, {ax, ay, az}, angle);
        vector<lld> ck = R({0, 0, 1}, {ax, ay, az}, angle);
        vector<vector<lld>> rotation_matrix = identity_matrix;
        rotation_matrix[0][0] = ci[0], rotation_matrix[1][0] = ci[1], rotation_matrix[2][0] = ci[2];
        rotation_matrix[0][1] = cj[0], rotation_matrix[1][1] = cj[1], rotation_matrix[2][1] = cj[2];
        rotation_matrix[0][2] = ck[0], rotation_matrix[1][2] = ck[1], rotation_matrix[2][2] = ck[2];
        v_stack[v_stack.size() - 1] = multiply(v_stack.back(), rotation_matrix);
    }
    // transforming a point
    vector<lld> transform(lld px, lld py, lld pz)
    {
        return transform_a_point(px, py, pz, v_stack.back());
    }
};

class view_transformation
{
private:
    vector<vector<lld>> transformation_matrix{4, vector<lld>(4, 0)};

public:
    // calculating the view transformation matrix
    view_transformation()
    {
        // V = RT
        // calculating the translation matrix
        vector<vector<lld>> translastion_matrix(4, vector<lld>(4, 0));
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                if (i == j)
                {
                    translastion_matrix[i][j] = 1;
                }
            }
        }
        translastion_matrix[0][3] = -eyeX;
        translastion_matrix[1][3] = -eyeY;
        translastion_matrix[2][3] = -eyeZ;

        // calculating rotation matrix
        vector<vector<lld>> rotation_matrix(4, vector<lld>(4, 0));
        rotation_matrix[3][3] = 1;
        vector<lld> up = {upX, upY, upZ};
        vector<lld> l = {lookX - eyeX, lookY - eyeY, lookZ - eyeZ};
        lld l_length = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
        for (int i = 0; i < 3; i++)
        {
            l[i] /= l_length;
        }
        vector<lld> r = cross(l, up);
        lld r_length = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        for (int i = 0; i < 3; i++)
        {
            r[i] /= r_length;
        }
        vector<lld> u = cross(r, l);
        rotation_matrix[0][0] = r[0];
        rotation_matrix[0][1] = r[1];
        rotation_matrix[0][2] = r[2];
        rotation_matrix[1][0] = u[0];
        rotation_matrix[1][1] = u[1];
        rotation_matrix[1][2] = u[2];
        rotation_matrix[2][0] = -l[0];
        rotation_matrix[2][1] = -l[1];
        rotation_matrix[2][2] = -l[2];
        transformation_matrix = multiply(rotation_matrix, translastion_matrix);
    }
    // transforming a point
    vector<lld> transform(lld px, lld py, lld pz)
    {
        return transform_a_point(px, py, pz, transformation_matrix);
    }
};

class projection_transformation
{
private:
    vector<vector<lld>> transformation_matrix{4, vector<lld>(4, 0)};

public:
    // calculating the projection matrix
    projection_transformation()
    {
        lld fovX = fovY * aspectRatio;
        lld t = near * tan((fovY / 2) * M_PI / 180L);
        lld r = near * tan((fovX / 2) * M_PI / 180L);
        transformation_matrix[0][0] = near / r;
        transformation_matrix[1][1] = near / t;
        transformation_matrix[2][2] = -(far + near) / (far - near);
        transformation_matrix[3][2] = -1;
        transformation_matrix[2][3] = -(2 * far * near) / (far - near);
    }
    // transforming a point
    vector<lld> transform(lld px, lld py, lld pz)
    {
        return transform_a_point(px, py, pz, transformation_matrix);
    }
};

void stage_1()
{
    ifstream fin("scene.txt");
    ofstream fout("stage1.txt");
    fout << fixed << setprecision(7);
    // taking input
    fin >> eyeX >> eyeY >> eyeZ;
    fin >> lookX >> lookY >> lookZ;
    fin >> upX >> upY >> upZ;
    fin >> fovY >> aspectRatio >> near >> far;
    transformation_stack ts;
    string cmd;
    while (fin >> cmd)
    {
        if (cmd == "triangle")
        {
            lld x1, y1, z1, x2, y2, z2, x3, y3, z3;
            fin >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
            vector<lld> p1 = ts.transform(x1, y1, z1);
            vector<lld> p2 = ts.transform(x2, y2, z2);
            vector<lld> p3 = ts.transform(x3, y3, z3);
            // outputing the transformed triangle
            fout << p1[0] << gap << p1[1] << gap << p1[2] << endl;
            fout << p2[0] << gap << p2[1] << gap << p2[2] << endl;
            fout << p3[0] << gap << p3[1] << gap << p3[2] << endl;
            fout << endl;
        }
        else if (cmd == "translate")
        {
            lld tx, ty, tz;
            fin >> tx >> ty >> tz;
            ts.translate(tx, ty, tz);
        }
        else if (cmd == "scale")
        {
            lld sx, sy, sz;
            fin >> sx >> sy >> sz;
            ts.scale(sx, sy, sz);
        }
        else if (cmd == "rotate")
        {
            lld angle, ax, ay, az;
            fin >> angle >> ax >> ay >> az;
            ts.rotate(angle, ax, ay, az);
        }
        else if (cmd == "push")
        {
            ts.push();
        }
        else if (cmd == "pop")
        {
            ts.pop();
        }
        else if (cmd == "end")
        {
            break;
        }
    }
    fin.close();
    fout.close();
}

void stage_2()
{
    ifstream fin("stage1.txt");
    ofstream fout("stage2.txt");
    fout << fixed << setprecision(7);
    view_transformation v;
    lld x1, y1, z1, x2, y2, z2, x3, y3, z3;
    while (fin >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3)
    {
        vector<lld> p1 = v.transform(x1, y1, z1);
        vector<lld> p2 = v.transform(x2, y2, z2);
        vector<lld> p3 = v.transform(x3, y3, z3);
        // outputing the transformed triangle
        fout << p1[0] << gap << p1[1] << gap << p1[2] << endl;
        fout << p2[0] << gap << p2[1] << gap << p2[2] << endl;
        fout << p3[0] << gap << p3[1] << gap << p3[2] << endl;
        fout << endl;
    }
    fin.close();
    fout.close();
}

void stage_3()
{
    ifstream fin("stage2.txt");
    ofstream fout("stage3.txt");
    fout << fixed << setprecision(7);
    projection_transformation pt;
    lld x1, y1, z1, x2, y2, z2, x3, y3, z3;
    while (fin >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3)
    {
        vector<lld> p1 = pt.transform(x1, y1, z1);
        vector<lld> p2 = pt.transform(x2, y2, z2);
        vector<lld> p3 = pt.transform(x3, y3, z3);
         // outputing the transformed triangle
        fout << p1[0] << gap << p1[1] << gap << p1[2] << endl;
        fout << p2[0] << gap << p2[1] << gap << p2[2] << endl;
        fout << p3[0] << gap << p3[1] << gap << p3[2] << endl;
        fout << endl;
    }
    fin.close();
    fout.close();
}

void z_buffer()
{
    ifstream con("config.txt");
    ifstream fin("stage3.txt");
    ofstream fout("z_buffer.txt");
    // taking input
    con >> Screen_Width >> Screen_Height;
    con >> left_limit_of_X;
    con >> bottom_limit_of_Y;
    con >> z_min >> z_max;
    right_limit_of_X = -left_limit_of_X, top_limit_of_Y = -bottom_limit_of_Y;
    // to map x and y coordinate pixel wise
    lld dx = (right_limit_of_X - left_limit_of_X) / Screen_Width;
    lld dy = (top_limit_of_Y - bottom_limit_of_Y) / Screen_Height;
    lld Top_Y = top_limit_of_Y - dy / 2;
    lld Left_X = left_limit_of_X + dx / 2;
    lld x1, y1, z1, x2, y2, z2, x3, y3, z3;
    // this matrix is initialized with infinite because z_value() function returns infinite for an invalid point
    vector<vector<lld>> z_buff(Screen_Width, vector<lld>(Screen_Height, infinite));
    // the image
    bitmap_image image(Screen_Width, Screen_Height);
    while (fin >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3)
    {
        vector<vector<lld>> triangle(3, vector<lld>(3, 0));
        triangle[0] = {x1, y1, z1};
        triangle[1] = {x2, y2, z2};
        triangle[2] = {x3, y3, z3};
        // for color
        int R = _random(), G = _random(), B = _random();

        // the (0,0) is the top left corner
        // the (Screen_Width - 1, Screen_Height - 1) is the bottom right corner

        // calculating the boundary box for the triangle and necessary clipping
        lld low_x = min({triangle[0][0], triangle[1][0], triangle[2][0]});
        low_x = max(left_limit_of_X, low_x);
        lld high_x = max({triangle[0][0], triangle[1][0], triangle[2][0]});
        high_x = min(right_limit_of_X, high_x);
        lld low_y = max({triangle[0][1], triangle[1][1], triangle[2][1]});
        low_y = min(top_limit_of_Y, low_y);
        lld high_y = min({triangle[0][1], triangle[1][1], triangle[2][1]});
        high_y = max(bottom_limit_of_Y, high_y);

        // calculating the range of row no and column no need to be checked
        int low_col = max(0, (int)floor((low_x - Left_X) / dx));
        int high_col = min(Screen_Width - 1, (int)ceil((high_x - Left_X) / dx));
        int low_row = max(0, (int)floor((Top_Y - low_y) / dy));
        int high_row = min(Screen_Height - 1, (int)ceil((Top_Y - high_y) / dy));

        // iterating
        for (int i = low_row; i <= high_row; i++)
        {
            for (int j = low_col; j <= high_col; j++)
            {
                // x, y value of a perticular pixel
                lld x = Left_X + j * dx;
                lld y = Top_Y - i * dy;
                // calculating the z value of the point
                lld z_val = z_value(triangle, {x, y});
                // checking the z value for range
                if (z_val >= z_min && z_val <= z_max)
                {
                    lld diff = z_buff[i][j] - z_val;
                    if (diff > 1e-8)
                    {
                        z_buff[i][j] = z_val;
                        image.set_pixel(j, i, R, G, B);
                    }
                }
            }
        }
    }
    // outputing z buffer value
    for (int i = 0; i < Screen_Width; i++)
    {
        for (int j = 0; j < Screen_Height; j++)
        {
            if (z_buff[i][j] != infinite)
                fout << z_buff[i][j] << '\t';
        }
        fout << endl;
    }
    image.save_image("out.bmp");
    con.close();
    fin.close();
    fout.close();
}

int main()
{
    stage_1();
    stage_2();
    stage_3();
    z_buffer();
    return 0;
}