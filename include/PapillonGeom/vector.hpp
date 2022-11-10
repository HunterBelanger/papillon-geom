/*
 * Papillon Geometry Library
 * Copyright 2022, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Papillon Geometry Library (PapillonGeom).
 *
 * PapillonGeom is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PapillonGeom is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PapillonGeom. If not, see <https://www.gnu.org/licenses/>.
 *
 * */
#ifndef PAPILLON_GEOM_VECTOR_H
#define PAPILLON_GEOM_VECTOR_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonGeom/pgeom_exception.hpp>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <ostream>

namespace pgeom {

class Direction;

/**
 * @brief Class which may represent any vector quantity such as position.
 */
class Vector {
 public:
  /**
   * @brief Construct a new Vector object with components (1,0,0).
   */
  Vector() : vals_{0., 0., 0.} {}

  /**
   * @brief Construct a new Vector object.
   * @param x  x component
   * @param y  y component
   * @param z  z component
   */
  Vector(double x, double y, double z) : vals_{x, y, z} {}

  /**
   * @brief Return a reference to desired component.
   * @param i Index of component, starting at 0.
   */
  double& operator[](std::size_t i) {
    switch (i) {
      case 0:
        return vals_[0];
        break;
      case 1:
        return vals_[1];
        break;
      case 2:
        return vals_[2];
        break;
      default:
        throw PGeomException("Index " + std::to_string(i) + "out of range.");
        return vals_[0];
        break;
    }
  }

  /**
   * @brief Return a const reference to desired component.
   * @param i  Index of component, starting at 0.
   */
  const double& operator[](std::size_t i) const {
    switch (i) {
      case 0:
        return vals_[0];
        break;
      case 1:
        return vals_[1];
        break;
      case 2:
        return vals_[2];
        break;
      default:
        throw PGeomException("Index " + std::to_string(i) + "out of range.");
        return vals_[0];
        break;
    }
  }

  /**
   * @brief Returns number of components (always 3).
   */
  constexpr std::size_t size() const { return 3; }

  /**
   * @brief Returns a pointer to the three components.
   */
  constexpr double* data() { return &vals_[0]; }

  /**
   * @brief Returns a const pointer to the three components.
   */
  constexpr const double* data() const { return data(); }

  /**
   * @brief Returns a reference to the x component.
   */
  constexpr double& x() { return vals_[0]; }

  /**
   * @brief Returns a const reference to the x component.
   */
  constexpr const double& x() const { return vals_[0]; }

  /**
   * @brief Returns a reference to the y component.
   */
  constexpr double& y() { return vals_[1]; }

  /**
   * @brief Returns a const reference to the y component.
   */
  constexpr const double& y() const { return vals_[1]; }

  /**
   * @brief Returns a reference to the z component.
   */
  constexpr double& z() { return vals_[2]; }

  /**
   * @brief Returns a const reference to the z component.
   */
  constexpr const double& z() const { return vals_[2]; }

  /**
   * @brief Returns the length or magnitude of the Vector.
   */
  double norm() const {
    return std::sqrt(vals_[0] * vals_[0] + vals_[1] * vals_[1] + vals_[2] * vals_[2]);
  }

  /**
   * @brief Returns the dot product with a second Vector.
   * @param v Vector with which to take the dot product.
   */
  double dot(const Vector& v) const {
    return vals_[0] * v.vals_[0] + vals_[1] * v.vals_[1] + vals_[2] * v.vals_[2];
  }

  /**
   * @brief Returns the dot product with a Direction.
   * @param d Direction with which to take the dot product.
   */
  double dot(const Direction& d) const;

  /**
   * @brief Returns the cross product with a second Vector.
   * @param v Vector with which to take the cross product.
   */
  Vector cross(const Vector& v) const {
    Vector out;
    out.x() = vals_[1] * v.vals_[2] - vals_[2] * v.vals_[1];
    out.y() = vals_[2] * v.vals_[0] - vals_[0] * v.vals_[2];
    out.z() = vals_[0] * v.vals_[1] - vals_[1] * v.vals_[0];
    return out;
  }

  /**
   * @brief Returns the cross product with a direction.
   * @param d Direction with which to take the cross product.
   */
  Vector cross(const Direction& d) const;

  /**
   * @brief Returns a copy of the original Vector.
   */
  Vector operator+() const { return *this; }

  /**
   * @brief Returns a Vector pointing in the oposite direction of this.
   */
  Vector operator-() const {
    return Vector(-this->x(), -this->y(), -this->z());
  }

  /**
   * @brief Returns the sum of two Vectors.
   * @param v Second vector in sum.
   */
  Vector operator+(const Vector& v) const {
    return Vector(vals_[0] + v.vals_[0], vals_[1] + v.vals_[1], vals_[2] + v.vals_[2]);
  }

  /**
   * @brief Returns difference of two Vectors.
   * @param v Vector to subtract from this.
   */
  Vector operator-(const Vector& v) const {
    return Vector(vals_[0] - v.vals_[0], vals_[1] - v.vals_[1], vals_[2] - v.vals_[2]);
  }

  /**
   * @brief Returns the dot product with a second Vector.
   * @param v Vector with which to take the dot product.
   */
  double operator*(const Vector& v) const {
    return dot(v); 
  }

  /**
   * @brief Returns a Vector which is this scaled by c.
   * @param c Scaling factor.
   */
  Vector operator*(const double& c) const {
    Vector out = *this;
    out.x() *= c;
    out.y() *= c;
    out.z() *= c;
    return out;
  }

  /**
   * @brief Returns a Vector which is this scaled by 1/c.
   * @param c Inverse of the scaling factor.
   */
  Vector operator/(const double& c) const {
    Vector out = *this;
    out.x() /= c;
    out.y() /= c;
    out.z() /= c;
    return out;
  }

  /**
   * @brief Addition assigment operator.
   * @param v Vector to add to this.
   */
  Vector& operator+=(const Vector& v) {
    for (std::size_t i = 0; i < 3; i++) {
      vals_[i] += v.vals_[i];
    }
    return *this;
  }

  /**
   * @brief Subtraction assigment operator.
   * @param v Vector to subtract from this.
   */
  Vector& operator-=(const Vector& v) {
    for (std::size_t i = 0; i < 3; i++) {
      vals_[i] -= v.vals_[i];
    }
    return *this;
  }

  /**
   * @brief Scalar multiplication assignment operator.
   * @param c Scalar factor.
   */
  Vector& operator*=(const double& c) {
    for (std::size_t i = 0; i < 3; i++) vals_[i] *= c;
    return *this;
  }

  /**
   * @brief Scalar divission assignment operator.
   * @param c Scalar factor.
   */
  Vector& operator/=(const double& c) {
    for (std::size_t i = 0; i < 3; i++) vals_[i] /= c;
    return *this;
  }

 protected:
  double vals_[3];
};

/**
 * @brief Returns a Vector which is this scaled by c.
 * @param v Vector to scale.
 * @param c Scaling factor.
 */
inline Vector operator*(const double& c, const Vector& v) {
  Vector out = v;
  out *= c;
  return out;
}

/**
 * @brief Writes a Vector to a stream.
 * @param strm Reference to std::ostream where data will be written.
 * @param v Vector to write.
 */
inline std::ostream& operator<<(std::ostream& strm, const Vector& v) {
  strm << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
  return strm;
}

/**
 * @brief Class which is a vector, representing the direction of travel in
 * 3D space. This vector is always guarenteed to be normalized to a magnitude
 * of one.
 */
class Direction : protected Vector {
 public:
  /**
   * @brief Construct a new Direction object with components (1,0,0).
   */
  Direction() : Vector(0., 0., 1.) {}

  /**
   * @brief Construct a new Direction object.
   * @param x x component
   * @param y y component
   * @param z z component
   */
  Direction(double x, double y, double z) : Vector(x, y, z) {
    double mag = this->norm();
    vals_[0] /= mag;
    vals_[1] /= mag;
    vals_[2] /= mag;
  }

  /**
   * @brief Construct a new Direction object.
   * @param mu
   * @param phi
   */
  Direction(double mu, double phi);

  /**
   * @brief Return a const reference to desired component.
   * @param i  Index of component, starting at 0.
   */
  const double& operator[](std::size_t i) const {
    return Vector::operator[](i); 
  }

  /**
   * @brief Returns a const pointer to the three components.
   */
  constexpr const double* data() const {
    return Vector::data(); 
  }

  /**
   * @brief Returns a const reference to the x component.
   */
  constexpr const double& x() const { return Vector::x(); }

  /**
   * @brief Returns a const reference to the y component.
   */
  constexpr const double& y() const { return Vector::y(); }

  /**
   * @brief Returns a const reference to the z component.
   */
  constexpr const double& z() const { return Vector::z(); }

  /**
   * @brief Returns the magnitude of the direction (in theory always 1).
   */
  double norm() const { return Vector::norm(); }

  /**
   * @brief Returns the dot product with a Vector.
   * @param v Vector with which to take the dot product.
   */
  double dot(const Vector& v) const {
    return v.x() * this->x() + v.y() * this->y() + v.z() * this->z();
  }

  /**
   * @brief Returns the dot product with a second Direction.
   * @param d Direction with which to take the dot product.
   */
  double dot(const Direction& d) const {
    return d.x() * this->x() + d.y() * this->y() + d.z() * this->z();
  }

  /**
   * @brief Returns the cross product with a Vector.
   * @param v Vector with which to take the cross product.
   */
  Vector cross(const Vector& v) const {
    Vector out;
    out.x() = vals_[1] * v[2] - vals_[2] * v[1];
    out.y() = vals_[2] * v[0] - vals_[0] * v[2];
    out.z() = vals_[0] * v[1] - vals_[1] * v[0];
    return out;
  }

  /**
   * @brief Returns the cross product with a second Direction.
   * @param d Direction with which to take the cross product.
   */
  Vector cross(const Direction& v) const {
    Vector out;
    out.x() = vals_[1] * v[2] - vals_[2] * v[1];
    out.y() = vals_[2] * v[0] - vals_[0] * v[2];
    out.z() = vals_[0] * v[1] - vals_[1] * v[0];
    return out;
  }

  /**
   * @brief Returns a copy of the original Direction
   */
  Direction operator+() const { return *this; }

  /**
   * @brief Returns a Direction pointing in the oposite direction of this.
   */
  Direction operator-() const {
    Direction out;
    out.vals_[0] = -vals_[0];
    out.vals_[1] = -vals_[1];
    out.vals_[2] = -vals_[2];
    return out;
  }

  /**
   * @brief Returns a Vector which is the sum of this and a Vector.
   * @param v Vector with which to take the sum.
   */
  Vector operator+(const Vector& v) const { return Vector::operator+(v); }

  /**
   * @brief Returns a Vector which is the sum of this and another Direction.
   * @param d Direction with which to take the sum.
   */
  Vector operator+(const Direction& d) const{ return Vector::operator+(d); }

  /**
   * @brief Returns a Vector which is the difference of this and a Vector.
   * @param v Vector to subtract from this.
   */
  Vector operator-(const Vector& v) const { return Vector::operator-(v); }

  /**
   * @brief Returns a Vector which is the difference of this and another
   * Direction.
   * @param d Direction to subtract from this.
   */
  Vector operator-(const Direction& d) const { return Vector::operator-(d); }

  /**
   * @brief Returns the dot product with a Vector.
   * @param v Vector with which to take the dot product.
   */
  double operator*(const Vector& v) const { return this->dot(v); }

  /**
   * @brief Returns the dot product with a second Direction.
   * @param d Direction with which to take the dot product.
   */
  double operator*(const Direction& d) const { return this->dot(d); }

  /**
   * @brief Returns a Vector which is this scaled by c.
   * @param c Scaling factor.
   */
  Vector operator*(const double& c) const { return Vector::operator*(c); }

  /**
   * @brief Returns a Vector which is this scaled by 1/c.
   * @param c Inverse of the scaling factor.
   */
  Vector operator/(const double& c) const { return Vector::operator/(c); }

};

/**
 * @brief Returns a Vector which is a Direction scaled by c.
 * @param c Scaling factor.
 * @param d Direction to scale.
 */
inline Vector operator*(const double& c, const Direction& d) {
  Vector out;
  out.x() = d.x() * c;
  out.y() = d.y() * c;
  out.z() = d.z() * c;
  return out;
}

/**
 * @brief Returns the dot product of a Vector with a Direction.
 * @param v Vector with which to take the dot product.
 * @param d Direction with which to take the dot product.
 */
inline double operator*(const Vector& v, const Direction& d) {
  return d.dot(v);
}

/**
 * @brief Returns a Vector which is the sum of a Vector and a direction.
 * @param v Vector with which to take the sum.
 * @param d Direction to add to v
 */
inline Vector operator+(const Vector& v, const Direction& d) { return d + v; }

/**
 * @brief Returns a Vector which is the difference of a Vector and a direction.
 * @param v Vector with which to take the difference.
 * @param d Direction to subtract from v
 */
inline Vector operator-(const Vector& v, const Direction& d) {
  return -(d - v);
}

/**
 * @brief Writes a Direction to a stream.
 * @param strm Reference to std::ostream where data will be written.
 * @param d Direction to write.
 */
inline std::ostream& operator<<(std::ostream& strm, const Direction& d) {
  strm << "<" << d.x() << ", " << d.y() << ", " << d.z() << ">";
  return strm;
}

inline double Vector::dot(const Direction& d) const {
  return vals_[0] * d[0] + vals_[1] * d[1] + vals_[2] * d[2];
}

inline Vector Vector::cross(const Direction& d) const {
  Vector out;
  out.x() = vals_[1] * d[2] - vals_[2] * d[1];
  out.y() = vals_[2] * d[0] - vals_[0] * d[2];
  out.z() = vals_[0] * d[1] - vals_[1] * d[0];
  return out;
}

using Position = Vector;

}  // namespace pgeom

#endif
