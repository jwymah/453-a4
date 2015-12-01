//---------------------------------------------------------------------------
//
// CPSC453 -- Introduction to Computer Graphics
// Assignment 2
//
// Classes and functions for manipulating points, vectors, matrices, 
// and colours.  You probably won't need to modify anything in these
// two files.
//
// Adapted from CS488 A2 University of Waterloo Computer Graphics Lab / 2003
//
//---------------------------------------------------------------------------

#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class Vector3D;

class Point2D
{
public:
  Point2D()
  {
    v_[0] = 0.0;
    v_[1] = 0.0;
  }
  Point2D(double x, double y)
  { 
    v_[0] = x;
    v_[1] = y;
  }
  Point2D(const Point2D& other)
  {
    v_[0] = other.v_[0];
    v_[1] = other.v_[1];
  }

  Point2D& operator =(const Point2D& other)
  {
    v_[0] = other.v_[0];
    v_[1] = other.v_[1];
    return *this;
  }

  double& operator[](size_t idx) 
  {
    return v_[ idx ];
  }
  double operator[](size_t idx) const 
  {
    return v_[ idx ];
  }

private:
  double v_[2];
};

class Point3D
{
public:
  Point3D()
  {
    v_[0] = 0.0;
    v_[1] = 0.0;
    v_[2] = 0.0;
  }
  Point3D(double x, double y, double z)
  {
    v_[0] = x;
    v_[1] = y;
    v_[2] = z;
  }
  Point3D(const Point3D& other)
  {
    v_[0] = other.v_[0];
    v_[1] = other.v_[1];
    v_[2] = other.v_[2];
  }

  Point3D& operator =(const Point3D& other)
  {
    v_[0] = other.v_[0];
    v_[1] = other.v_[1];
    v_[2] = other.v_[2];
    return *this;
  }

  double& operator[](size_t idx)
  {
    return v_[ idx ];
  }
  double operator[](size_t idx) const
  {
    return v_[ idx ];
  }

private:
  double v_[3];
};

class Point4D
{
public:
  Point4D()
  {
    v_[0] = 0.0;
    v_[1] = 0.0;
    v_[2] = 0.0;
    v_[3] = 0.0;
  }
  Point4D(double x, double y, double z, double w)
  {
    v_[0] = x;
    v_[1] = y;
    v_[2] = z;
    v_[3] = w;
  }
  Point4D(const Point4D& other)
  {
    v_[0] = other.v_[0];
    v_[1] = other.v_[1];
    v_[2] = other.v_[2];
    v_[3] = other.v_[3];
  }

  Point4D& operator =(const Point4D& other)
  {
    v_[0] = other.v_[0];
    v_[1] = other.v_[1];
    v_[2] = other.v_[2];
    v_[3] = other.v_[3];
    return *this;
  }

  double& operator[](size_t idx)
  {
    return v_[ idx ];
  }
  double operator[](size_t idx) const
  {
    return v_[ idx ];
  }

private:
  double v_[4];
};

class Vector3D
{
public:
  Vector3D()
  {
    v_[0] = 0.0;
    v_[1] = 0.0;
    v_[2] = 0.0;
  }
  Vector3D(double x, double y, double z)
  {
    v_[0] = x;
    v_[1] = y;
    v_[2] = z;
  }
  Vector3D(const Vector3D& other)
  {
    v_[0] = other.v_[0];
    v_[1] = other.v_[1];
    v_[2] = other.v_[2];
  }

  Vector3D& operator =(const Vector3D& other)
  {
    v_[0] = other.v_[0];
    v_[1] = other.v_[1];
    v_[2] = other.v_[2];
    return *this;
  }

  Vector3D operator *(const Vector3D& other)
  {
      Vector3D newVec = Vector3D(
    v_[0] * other.v_[0],
    v_[1] * other.v_[1],
    v_[2] * other.v_[2]);
      return newVec;
  }

  double& operator[](size_t idx)
  {
    return v_[ idx ];
  }
  double operator[](size_t idx) const
  {
    return v_[ idx ];
  }

  double dot(const Vector3D& other) const
  {
    return v_[0]*other.v_[0] + v_[1]*other.v_[1] + v_[2]*other.v_[2];
  }

  double length2() const
  {
    return v_[0]*v_[0] + v_[1]*v_[1] + v_[2]*v_[2];
  }
  double length() const
  {
    return sqrt(length2());
  }

  double normalize();

  Vector3D cross(const Vector3D& other) const
  {
    return Vector3D(
                    v_[1]*other[2] - v_[2]*other[1],
                    v_[2]*other[0] - v_[0]*other[2],
                    v_[0]*other[1] - v_[1]*other[0]);
  }

private:
  double v_[3];
};

inline Vector3D operator *(double s, const Vector3D& v)
{
  return Vector3D(s*v[0], s*v[1], s*v[2]);
}

inline Vector3D operator +(const Vector3D& a, const Vector3D& b)
{
  return Vector3D(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}

inline Point3D operator +(const Point3D& a, const Vector3D& b)
{
  return Point3D(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}

inline Vector3D operator -(const Point3D& a, const Point3D& b)
{
  return Vector3D(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}

inline Vector3D operator -(const Vector3D& a, const Vector3D& b)
{
  return Vector3D(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}

inline Vector3D operator -(const Vector3D& a)
{
  return Vector3D(-a[0], -a[1], -a[2]);
}

inline Point3D operator -(const Point3D& a, const Vector3D& b)
{
  return Point3D(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}

inline Vector3D cross(const Vector3D& a, const Vector3D& b) 
{
  return a.cross(b);
}

inline std::ostream& operator <<(std::ostream& os, const Point2D& p)
{
  return os << "p<" << p[0] << "," << p[1] << ">";
}

inline std::ostream& operator <<(std::ostream& os, const Point3D& p)
{
  return os << "p<" << p[0] << "," << p[1] << "," << p[2] << ">";
}

inline std::ostream& operator <<(std::ostream& os, const Vector3D& v)
{
  return os << "v<" << v[0] << "," << v[1] << "," << v[2] << ">";
}

class Matrix4x4;

class Vector4D
{
public:
  Vector4D()
  {
    v_[0] = 0.0;
    v_[1] = 0.0;
    v_[2] = 0.0;
    v_[3] = 0.0;
  }
  Vector4D(double x, double y, double z, double w)
  { 
    v_[0] = x;
    v_[1] = y;
    v_[2] = z;
    v_[3] = w;
  }

  Vector4D(const Vector4D& other)
  {
    v_[0] = other.v_[0];
    v_[1] = other.v_[1];
    v_[2] = other.v_[2];
    v_[3] = other.v_[3];
  }

  Vector4D& operator =(const Vector4D& other)
  {
    v_[0] = other.v_[0];
    v_[1] = other.v_[1];
    v_[2] = other.v_[2];
    v_[3] = other.v_[3];
    return *this;
  }

  double& operator[](size_t idx) 
  {
    return v_[ idx ];
  }
  double operator[](size_t idx) const 
  {
    return v_[ idx ];
  }

private:
  double v_[4];
};

class Matrix4x4
{
public:
  Matrix4x4()
  {
    // Construct an identity matrix
    std::fill(v_, v_+16, 0.0);
    v_[0] = 1.0;
    v_[5] = 1.0;
    v_[10] = 1.0;
    v_[15] = 1.0;
  }
  Matrix4x4(const Matrix4x4& other)
  {
    std::copy(other.v_, other.v_+16, v_);
  }
  Matrix4x4(const Vector4D row1, const Vector4D row2, const Vector4D row3, 
             const Vector4D row4)
  {
    v_[0] = row1[0]; 
    v_[1] = row1[1]; 
    v_[2] = row1[2]; 
    v_[3] = row1[3]; 

    v_[4] = row2[0]; 
    v_[5] = row2[1]; 
    v_[6] = row2[2]; 
    v_[7] = row2[3]; 

    v_[8] = row3[0]; 
    v_[9] = row3[1]; 
    v_[10] = row3[2]; 
    v_[11] = row3[3]; 

    v_[12] = row4[0]; 
    v_[13] = row4[1]; 
    v_[14] = row4[2]; 
    v_[15] = row4[3]; 
  }
  Matrix4x4(double *vals)
  {
    std::copy(vals, vals + 16, (double*)v_);
  }

  Matrix4x4& operator=(const Matrix4x4& other)
  {
    std::copy(other.v_, other.v_+16, v_);
    return *this;
  }

  Vector4D getRow(size_t row) const
  {
    return Vector4D(v_[4*row], v_[4*row+1], v_[4*row+2], v_[4*row+3]);
  }
  double *getRow(size_t row) 
  {
    return (double*)v_ + 4*row;
  }

  Vector4D getColumn(size_t col) const
  {
    return Vector4D(v_[col], v_[4+col], v_[8+col], v_[12+col]);
  }

  Vector4D operator[](size_t row) const
  {
    return getRow(row);
  }
  double *operator[](size_t row) 
  {
    return getRow(row);
  }

  Matrix4x4 transpose() const
  {
    return Matrix4x4(getColumn(0), getColumn(1), 
                      getColumn(2), getColumn(3));
  }
  Matrix4x4 invert() const;

  const double *begin() const
  {
    return (double*)v_;
  }
  const double *end() const
  {
    return begin() + 16;
  }
		
private:
  double v_[16];
};

inline Matrix4x4 operator *(const Matrix4x4& a, const Matrix4x4& b)
{
  Matrix4x4 ret;

  for(size_t i = 0; i < 4; ++i) {
    Vector4D row = a.getRow(i);
		
    for(size_t j = 0; j < 4; ++j) {
      ret[i][j] = row[0] * b[0][j] + row[1] * b[1][j] + 
        row[2] * b[2][j] + row[3] * b[3][j];
    }
  }

  return ret;
}

inline Vector3D operator *(const Matrix4x4& M, const Vector3D& v)
{
  return Vector3D(
                  v[0] * M[0][0] + v[1] * M[0][1] + v[2] * M[0][2],
                  v[0] * M[1][0] + v[1] * M[1][1] + v[2] * M[1][2],
                  v[0] * M[2][0] + v[1] * M[2][1] + v[2] * M[2][2]);
}

inline Point3D operator *(const Matrix4x4& M, const Point3D& p)
{
  return Point3D(
                 p[0] * M[0][0] + p[1] * M[0][1] + p[2] * M[0][2] + M[0][3],
                 p[0] * M[1][0] + p[1] * M[1][1] + p[2] * M[1][2] + M[1][3],
                 p[0] * M[2][0] + p[1] * M[2][1] + p[2] * M[2][2] + M[2][3]);
}

inline Vector3D transNorm(const Matrix4x4& M, const Vector3D& n)
{
  return Vector3D(
                  n[0] * M[0][0] + n[1] * M[1][0] + n[2] * M[2][0],
                  n[0] * M[0][1] + n[1] * M[1][1] + n[2] * M[2][1],
                  n[0] * M[0][2] + n[1] * M[1][2] + n[2] * M[2][2]);
}

inline std::ostream& operator <<(std::ostream& os, const Matrix4x4& M)
{
  return os << "[" << M[0][0] << " " << M[0][1] << " " 
            << M[0][2] << " " << M[0][3] << "]" << std::endl
            << "[" << M[1][0] << " " << M[1][1] << " " 
            << M[1][2] << " " << M[1][3] << "]" << std::endl
            << "[" << M[2][0] << " " << M[2][1] << " " 
            << M[2][2] << " " << M[2][3] << "]" << std::endl
            << "[" << M[3][0] << " " << M[3][1] << " " 
            << M[3][2] << " " << M[3][3] << "]";
}

class Colour
{
public:
  Colour(double r, double g, double b)
    : r_(r)
    , g_(g)
    , b_(b)
  {}
  Colour(double c)
    : r_(c)
    , g_(c)
    , b_(c)
  {}
  Colour(const Colour& other)
    : r_(other.r_)
    , g_(other.g_)
    , b_(other.b_)
  {}

  Colour& operator =(const Colour& other)
  {
    r_ = other.r_;
    g_ = other.g_;
    b_ = other.b_;
    return *this;
  }

  double R() const 
  { 
    return r_;
  }
  double G() const 
  { 
    return g_;
  }
  double B() const 
  { 
    return b_;
  }

private:
  double r_;
  double g_;
  double b_;
};

inline Colour operator *(double s, const Colour& a)
{
  return Colour(s*a.R(), s*a.G(), s*a.B());
}

inline Colour operator *(const Colour& a, const Colour& b)
{
  return Colour(a.R()*b.R(), a.G()*b.G(), a.B()*b.B());
}

inline Colour operator +(const Colour& a, const Colour& b)
{
  return Colour(a.R()+b.R(), a.G()+b.G(), a.B()+b.B());
}

inline std::ostream& operator <<(std::ostream& os, const Colour& c)
{
  return os << "c<" << c.R() << "," << c.G() << "," << c.B() << ">";
}

class Line3D
{
public:
    Line3D();
    ~Line3D();
    Line3D(const Point3D &p1, Point3D &p2);
    Line3D(const Line3D &other);
    Point3D &getP1();
    Point3D &getP2();
private:
    Point3D p1_, p2_;
};

class Triangle
{
public:
    Triangle();
    ~Triangle();

    Matrix4x4 getTransform() const;
    void resetTransform();
    void appendTransform(const Matrix4x4 &xform);

    std::vector<Line3D> getLines();

private:
    std::vector<Point3D> verts_;
    Matrix4x4 transform_;
};

class Gnomon
{
public:
    Matrix4x4 getTransform() const;
    void resetTransform();
    void appendTranslationTransform(const Matrix4x4 &xform);
    void appendRotationTransform(const Matrix4x4 &xform);

    std::vector<Line3D> getLines();
    Gnomon();
    ~Gnomon();

private:
    std::vector<Point3D> verts_;
    Matrix4x4 translationTransform_;
    Matrix4x4 rotationTransform_;
};

class Cube
{
public:
    Cube();
    ~Cube();

    Matrix4x4 getTransform() const; //TODO: possibly split this up to get individual xforms
    Gnomon gnomon;
    void resetTransform();
    void appendTranslationTransform(const Matrix4x4 &xform);
    void appendRotationTransform(const Matrix4x4 &xform);

    float scale_factor_x;
    float scale_factor_y;
    float scale_factor_z;

    void setScaleX(float scale);
    void setScaleY(float scale);
    void setScaleZ(float scale);

    Matrix4x4 getScaleTransform() const;

    std::vector<Line3D> getLines();

private:
    std::vector<Point3D> verts_;
    Matrix4x4 translationTransform_;
    Matrix4x4 rotationTransform_;
    int edgeLength;
};





class Shape
{
public:
    Vector3D surfaceColour, kd, ks;
    float transparency;
    bool intersect(Point3D rayOrigin, Vector3D rayDir)
    {
        float t0, t1;
        return intersect(rayOrigin, rayDir, t0, t1);
    }

    virtual bool intersect(Point3D rayOrigin, Vector3D rayDir, float &t0, float &t1) = 0;
    // assumes intersecti0on
    virtual Vector3D normalAt(Point3D p) = 0;
};

class Sphere: public Shape
{
public:
    Point3D center;
    float radius, radius2;

    Sphere(Point3D c, float r, Vector3D sc, Vector3D kdVector, Vector3D ksVector, float transp = 0);
    ~Sphere();

    bool intersect(Point3D rayOrigin, Vector3D rayDir, float &t0, float &t1)
    {
        float rayOriginX = rayOrigin[0];
        float rayOriginY = rayOrigin[1];
        float rayOriginZ = rayOrigin[2];
        float rayDirX = rayDir[0];
        float rayDirY = rayDir[1];
        float rayDirZ = rayDir[2];
        float centerX = center[0];
        float centerY = center[1];
        float centerZ = center[2];

        float a = pow(rayDirX, 2) + pow(rayDirY, 2) + pow(rayDirZ, 2); //don't need this because rayDir is normalized
        float b = 2 * (rayDirX * (rayOriginX - centerX)
                 + rayDirY * (rayOriginY - centerY)
                 + rayDirZ * (rayOriginZ - centerZ)
                 );
        float c = pow(rayOriginX - centerX, 2)
                + pow(rayOriginY - centerY, 2)
                + pow(rayOriginZ - centerZ, 2)
                - radius2;

        t0 = (-b - sqrt(pow(b, 2) - 4*c)) / 2;
        t1 = (-b + sqrt(pow(b, 2) - 4*c)) / 2;

        if (t0 > 0.0)
        {
            return true;
        }
        if (t1 > 0.0)
        {
            return true;
        }
        return false;

////       Geometric solution, not competely sure how it works
//        Vector3D l = center - rayOrigin;
//        float tca = l.dot(rayDir);
//        if (tca < 0)
//        {
//            return false;
//        }
//        float d2 = l.dot(l) - tca * tca;
//        if (d2 > radius2)
//        {
//            return false;
//        }
//        return true;
    }

    Vector3D normalAt(Point3D p)
    {
        Vector3D n = Vector3D(p[0] - center[0],
                        p[1] - center[1],
                        p[2] - center[2]);
        n.normalize();
        return n;
    }
};

#endif // ALGEBRA_H

