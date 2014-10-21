

#pragma once

#include <math.h>
#include <vector>
#include <string>

#include "lighter.h"


#ifdef _MSC_VER
#define FORCEINLINE __forceinline
#else
#define FORCEINLINE static inline __attribute__((__always_inline__))
#endif


#define LTR_WT_COLINFO  0 // generate collision info (* mesh)
#define LTR_WT_SAMPLES  1 // gather sampled information to compact position/normal/index arrays, allocate color arrays (* mesh)
#define LTR_WT_LMRENDER 2 // generate colors for samples by rendering lights to lightmaps (* mesh * light)
#define LTR_WT_FINALIZE 3 // push sample data back to lightmaps (* mesh)


template< class T > T TMIN( const T& a, const T& b ){ return a < b ? a : b; }
template< class T > T TMAX( const T& a, const T& b ){ return a > b ? a : b; }
template< class T > void TMEMSET( T* a, size_t c, const T& v )
{
	for( size_t i = 0; i < c; ++i )
		a[ i ] = v;
}


struct Vec2
{
	float x, y;
	
	static Vec2 Create( float x, float y ){ Vec2 v = { x, y }; return v; }
	
	Vec2 operator + ( const Vec2& o ) const { Vec2 v = { x + o.x, y + o.y }; return v; }
	Vec2 operator - ( const Vec2& o ) const { Vec2 v = { x - o.x, y - o.y }; return v; }
	Vec2 operator * ( const Vec2& o ) const { Vec2 v = { x * o.x, y * o.y }; return v; }
	Vec2 operator / ( const Vec2& o ) const { Vec2 v = { x / o.x, y / o.y }; return v; }
	
	Vec2 operator + ( float f ) const { Vec2 v = { x + f, y + f }; return v; }
	Vec2 operator - ( float f ) const { Vec2 v = { x - f, y - f }; return v; }
	Vec2 operator * ( float f ) const { Vec2 v = { x * f, y * f }; return v; }
	Vec2 operator / ( float f ) const { Vec2 v = { x / f, y / f }; return v; }
	
	Vec2& operator += ( const Vec2& o ){ x += o.x; y += o.y; return *this; }
	Vec2& operator -= ( const Vec2& o ){ x -= o.x; y -= o.y; return *this; }
	Vec2& operator *= ( const Vec2& o ){ x *= o.x; y *= o.y; return *this; }
	Vec2& operator /= ( const Vec2& o ){ x /= o.x; y /= o.y; return *this; }
	
	Vec2& operator += ( float f ){ x += f; y += f; return *this; }
	Vec2& operator -= ( float f ){ x -= f; y -= f; return *this; }
	Vec2& operator *= ( float f ){ x *= f; y *= f; return *this; }
	Vec2& operator /= ( float f ){ x /= f; y /= f; return *this; }
};

FORCEINLINE float Vec2CrossProduct( const Vec2& v1, const Vec2& v2 )
{
	return ( v1.x * v2.y ) - ( v1.y * v2.x );
}

struct Vec3
{
	float x, y, z;
	
	static Vec3 Create( float x, float y, float z ){ Vec3 v = { x, y, z }; return v; }
	
	Vec3 operator + ( const Vec3& o ) const { Vec3 v = { x + o.x, y + o.y, z + o.z }; return v; }
	Vec3 operator - ( const Vec3& o ) const { Vec3 v = { x - o.x, y - o.y, z - o.z }; return v; }
	Vec3 operator * ( const Vec3& o ) const { Vec3 v = { x * o.x, y * o.y, z * o.z }; return v; }
	Vec3 operator / ( const Vec3& o ) const { Vec3 v = { x / o.x, y / o.y, z / o.z }; return v; }
	
	Vec3 operator + ( float f ) const { Vec3 v = { x + f, y + f, z + f }; return v; }
	Vec3 operator - ( float f ) const { Vec3 v = { x - f, y - f, z - f }; return v; }
	Vec3 operator * ( float f ) const { Vec3 v = { x * f, y * f, z * f }; return v; }
	Vec3 operator / ( float f ) const { Vec3 v = { x / f, y / f, z / f }; return v; }
	
	Vec3& operator += ( const Vec3& o ){ x += o.x; y += o.y; z += o.z; return *this; }
	Vec3& operator -= ( const Vec3& o ){ x -= o.x; y -= o.y; z -= o.z; return *this; }
	Vec3& operator *= ( const Vec3& o ){ x *= o.x; y *= o.y; z *= o.z; return *this; }
	Vec3& operator /= ( const Vec3& o ){ x /= o.x; y /= o.y; z /= o.z; return *this; }
	
	Vec3& operator += ( float f ){ x += f; y += f; z += f; return *this; }
	Vec3& operator -= ( float f ){ x -= f; y -= f; z -= f; return *this; }
	Vec3& operator *= ( float f ){ x *= f; y *= f; z *= f; return *this; }
	Vec3& operator /= ( float f ){ x /= f; y /= f; z /= f; return *this; }
};

struct Mat4
{
	union
	{
		float a[16];
		float m[4][4];
	};
};


typedef std::vector< u32 > U32Vector;
typedef std::vector< float > FloatVector;
typedef std::vector< Vec2 > Vec2Vector;
typedef std::vector< Vec3 > Vec3Vector;
typedef std::vector< Mat4 > Mat4Vector;
typedef std::vector< ltr_LightInfo > LightVector;
typedef std::vector< ltr_WorkOutput > WorkOutputVector;


void RasterizeTriangle2D( Vec3* image, i32 width, i32 height, const Vec2& p1, const Vec2& p2, const Vec2& p3, const Vec3& v1, const Vec3& v2, const Vec3& v3 );


