#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <math.h>

/// @brief 3 dimensional vector structure.
struct vector3d
{
    float x;
    float y;
    float z;

public:
    vector3d();
    vector3d(float in_x, float in_y, float in_z);
    vector3d(float inxyz);
    ~vector3d();

    /// @brief Get the x-axis component
    /// @return
    float X() { return x; }
    /// @brief Get the y-axis component
    /// @return
    float Y() { return y; }
    /// @brief Get the z-axis component
    /// @return
    float Z() { return z; }
    /// @brief Return the magnitude/length of this vector
    /// @return
    float Magnitude()
    {
        return sqrtf((x * x) + (y * y) + (z * z));
    }
    /// @brief  Gets the squared length of the vector.
    /// @return
    float SquareMagnitude()
    {
        return ((x * x) + (y * y) + (z * z));
    }
    /// @brief Get the normal of this vector
    /// @return
    vector3d Normal()
    {
        float length = Magnitude();
        if (length > 0)
        {
            return { x / length, y / length, z / length };
        }
        else
        {
            return vector3d();
        }
    }
    /// @brief Scale this vector components to unit sizes
    void Normalise()
    {
        float length = Magnitude();
        if (length > 0)
        {
            x /= length;
            y /= length;
            z /= length;
        }
    }
    /// @brief Normalise the rotation angles
    inline void NormaliseEuler()
    {
    }
    /// @brief Check if any of the vector components are out of range
    /// @return
    bool IsNaN()
    {
        return isnan(x) || isnan(y) || isnan(z);
    }
    /// @brief Check if any of the vector components are out of range
    /// @return
    bool IsInf()
    {
        return isinf(x) || isinf(y) || isinf(z);
    }
    /// @brief Zero all the components of the vector.
    void Clear()
    {
        x = 0.0f;
        y = 0.0f;
        z = 0.0f;
    }
    /// @brief Get the cross product of two vectors
    /// @param rhs
    /// @param lhs
    static inline vector3d CrossProduct(vector3d lhs, vector3d rhs)
    {
        return {
            ((lhs.y * rhs.z) - (lhs.z * rhs.y)),
            ((lhs.z * rhs.x) - (lhs.x * rhs.z)),
            ((lhs.x * rhs.y) - (lhs.y * rhs.x)) };
    }
    /// @brief Get the dot product of two vectors
    /// @param lhs
    /// @param rhs
    /// @return
    static inline float DotProduct(vector3d lhs, vector3d rhs)
    {
        return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
    }

    static vector3d Zero() { return vector3d(0, 0, 0); }
    /// @brief A vector with components (0,0,0);
    /// @return
    static vector3d Left() { return vector3d(1, 0, 0); }
    /// @brief A vector with components (-1,0,0);
    /// @return
    static vector3d Right() { return vector3d(-1, 0, 0); }
    /// @brief A vector with components (0,1,0);
    /// @return
    static vector3d Up() { return vector3d(0, 1, 0); }
    /// @brief A vector with components (0,-1,0);
    /// @return
    static vector3d Down() { return vector3d(0, -1, 0); }
    /// @brief A vector with components (0,0,1);
    /// @return
    static vector3d Backward() { return vector3d(0, 0, 1); }
    /// @brief A vector with components (0,0,-1);
    /// @return
    static vector3d Forward() { return vector3d(0, 0, -1); }
    /// @brief A vector with components (1,1,1);
    /// @return
    static vector3d One() { return vector3d(1, 1, 1); }
};

/// @brief Base constructor
inline vector3d::vector3d()
{
    x = 0;
    y = 0;
    z = 0;
}
/// @brief construct a vector with the specified components
/// @param in_x
/// @param in_y
/// @param in_z
inline vector3d::vector3d(float in_x, float in_y, float in_z)
{
    x = in_x;
    y = in_y;
    z = in_z;
}
/// @brief Construct a vector with the components set to the specified value
/// @param inxyz
inline vector3d::vector3d(float inxyz)
{
    x = inxyz;
    y = inxyz;
    z = inxyz;
}
/// @brief Destructor
inline vector3d::~vector3d() {}

/// @brief Multiply vector with float
/// @param lhs
/// @param rhs
/// @return
static inline vector3d operator*(vector3d lhs, float rhs)
{
    return vector3d(
        lhs.x * rhs,
        lhs.y * rhs,
        lhs.z * rhs);
}
/// @brief Multiply vector with float
/// @param lhs
/// @param rhs
/// @return
static inline vector3d operator*(float rhs, vector3d lhs)
{
    return vector3d(
        lhs.x * rhs,
        lhs.y * rhs,
        lhs.z * rhs);
}
/// @brief Add a float to each component of the specified vector
/// @param lhs
/// @param rhs
/// @return
static inline vector3d operator+(vector3d lhs, float rhs)
{
    return vector3d(
        lhs.x + rhs,
        lhs.y + rhs,
        lhs.z + rhs);
}
/// @brief Add a float to each component of the specified vector
/// @param lhs
/// @param rhs
/// @return
static inline vector3d operator-(vector3d lhs, float rhs)
{
    return vector3d(
        lhs.x - rhs,
        lhs.y - rhs,
        lhs.z - rhs);
}
/// @brief Divide vector by float
/// @param lhs
/// @param rhs
/// @return
static inline vector3d operator/(vector3d lhs, float rhs)
{
    return vector3d(
        lhs.x / rhs,
        lhs.y / rhs,
        lhs.z / rhs);
}
/// @brief Unary minus
/// @param value
/// @return
static inline vector3d operator-(vector3d value)
{
    return vector3d(-value.x, -value.y, -value.z);
}

/// @brief Add two vectors together
/// @param lhs
/// @param rhs
/// @return
static inline vector3d operator+(vector3d lhs, vector3d rhs)
{
    return vector3d(
        lhs.x + rhs.x,
        lhs.y + rhs.y,
        lhs.z + rhs.z);
}
/// @brief Get the difference of two vectors
/// @param lhs
/// @param rhs
/// @return
static inline vector3d operator-(vector3d lhs, vector3d rhs)
{
    return vector3d(
        lhs.x - rhs.x,
        lhs.y - rhs.y,
        lhs.z - rhs.z);
}

#endif