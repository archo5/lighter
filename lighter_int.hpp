

#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <string>

#include "lighter.h"


#ifdef _WIN32

#  define WIN32_LEAN_AND_MEAN
#  undef _WIN32_WINNT
#  define _WIN32_WINNT 0x0600
#  undef WINVER
#  define WINVER 0x0600
#  include <windows.h>

#  define ltrthread_sleep( ms ) Sleep( (DWORD) ms )

#  define ltrmutex_t CRITICAL_SECTION
#  define ltrthread_t HANDLE
#  define threadret_t DWORD __stdcall
#  define threadarg_t void*

#  define ltrthread_create( toT, func, data ) toT = CreateThread( NULL, 1024, func, data, 0, NULL )
#  define ltrthread_self() GetCurrentThread()
#  define ltrthread_join( T ) for(;;){ WaitForMultipleObjects( 1, &T, TRUE, INFINITE ); CloseHandle( T ); break; }
#  define ltrthread_equal( T1, T2 ) (T1 == T2)

#  define ltrmutex_init( M ) InitializeCriticalSection( &M )
#  define ltrmutex_destroy( M ) DeleteCriticalSection( &M )
#  define ltrmutex_lock( M ) EnterCriticalSection( &M )
#  define ltrmutex_unlock( M ) LeaveCriticalSection( &M )

static int ltrnumcpus()
{
	SYSTEM_INFO sysinfo;
	GetSystemInfo( &sysinfo );
	return sysinfo.dwNumberOfProcessors;
}


#else

#  include <unistd.h>
#  include <pthread.h>

static void ltrthread_sleep( uint32_t ms )
{
	if( ms >= 1000 )
	{
		sleep( ms / 1000 );
		ms %= 1000;
	}
	if( ms > 0 )
	{
		usleep( ms * 1000 );
	}
}

#  define ltrmutex_t pthread_mutex_t
#  define ltrthread_t pthread_t
#  define threadret_t void*
#  define threadarg_t void*

#  define ltrthread_create( toT, func, data ) pthread_create( &toT, NULL, func, data )
#  define ltrthread_self() pthread_self()
#  define ltrthread_join( T ) pthread_join( T, NULL )
#  define ltrthread_equal( T1, T2 ) pthread_equal( T1, T2 )

#  define ltrmutex_init( M ) pthread_mutex_init( &M, NULL )
#  define ltrmutex_destroy( M ) pthread_mutex_destroy( &M )
#  define ltrmutex_lock( M ) pthread_mutex_lock( &M )
#  define ltrmutex_unlock( M ) pthread_mutex_unlock( &M )

static int ltrnumcpus()
{
	return sysconf( _SC_NPROCESSORS_ONLN );
}


#endif


#ifndef FORCEINLINE
#ifdef _MSC_VER
#define FORCEINLINE __forceinline
#else
#define FORCEINLINE inline __attribute__((__always_inline__))
#endif
#endif


struct LTRMutex
{
	LTRMutex(){ ltrmutex_init( mutex ); }
	~LTRMutex(){ ltrmutex_destroy( mutex ); }
	void Lock(){ ltrmutex_lock( mutex ); }
	void Unlock(){ ltrmutex_unlock( mutex ); }
	void SleepCV( CONDITION_VARIABLE* cv ){ SleepConditionVariableCS( cv, &mutex, INFINITE ); }
	ltrmutex_t mutex;
};
struct LTRMutexLock
{
	LTRMutexLock( LTRMutex& m ) : mutex( m ){ m.Lock(); }
	~LTRMutexLock(){ mutex.Unlock(); }
	LTRMutex& mutex;
};

struct LTRWorker
{
	struct IO
	{
		void* shared;
		void* item;
		size_t i;
	};
	typedef void (*WorkProc) (IO*);
	
	LTRWorker() :
		m_shared( NULL ),
		m_items( NULL ),
		m_itemSize( 0 ),
		m_itemCount( 0 ),
		m_nextItem( 0 ),
		m_numDone( 0 ),
		m_workProc( NULL ),
		m_exit( false )
	{
		InitializeConditionVariable( &m_hasWork );
		InitializeConditionVariable( &m_hasDone );
	}
	~LTRWorker()
	{
		m_exit = true;
		WakeAllConditionVariable( &m_hasWork );
		for( size_t i = 0; i < m_threads.size(); ++i )
		{
			ltrthread_join( m_threads[ i ] );
		}
	}
	void Init( int nt = ltrnumcpus() )
	{
		if( m_threads.size() )
			return;
		m_threads.resize( nt );
		for( int i = 0; i < nt; ++i )
		{
			ltrthread_create( m_threads[ i ], threadproc, this );
		}
	}
	
	void DoWork( void* shared, void* items, size_t size, size_t count, WorkProc wp )
	{
		m_mutex.Lock();
		
		m_shared = shared;
		m_items = (char*) items;
		m_itemSize = size;
		m_itemCount = count;
		m_nextItem = 0;
		m_numDone = 0;
		m_workProc = wp;
		
		m_mutex.Unlock();
		WakeAllConditionVariable( &m_hasWork );
		m_mutex.Lock();
		
		while( m_numDone < m_itemCount )
		{
			m_mutex.SleepCV( &m_hasDone );
		}
		
		m_shared = NULL;
		m_items = NULL;
		m_itemSize = 0;
		m_itemCount = 0;
		m_nextItem = 0;
		m_numDone = 0;
		m_workProc = NULL;
		
		m_mutex.Unlock();
	}
	
	void IntProcess()
	{
		m_mutex.Lock();
		while( !m_exit )
		{
			if( m_nextItem < m_itemCount )
			{
				size_t myitem = m_nextItem++;
				m_mutex.Unlock();
				
				IO io = { m_shared, m_items + myitem * m_itemSize, myitem };
				m_workProc( &io );
				
				m_mutex.Lock();
				m_numDone++;
				if( m_numDone < m_itemCount )
					continue; // there may be more work
			}
			
			m_mutex.Unlock();
			WakeAllConditionVariable( &m_hasDone );
			m_mutex.Lock();
			
			m_mutex.SleepCV( &m_hasWork );
		}
		m_mutex.Unlock();
	}
	
	static threadret_t threadproc( threadarg_t arg )
	{
		LTRWorker* w = (LTRWorker*) arg;
		w->IntProcess();
		return 0;
	}
	
	void* m_shared;
	char* m_items;
	size_t m_itemSize;
	size_t m_itemCount;
	size_t m_nextItem;
	size_t m_numDone;
	WorkProc m_workProc;
	
	std::vector< ltrthread_t > m_threads;
	LTRMutex m_mutex;
	volatile bool m_exit;
	CONDITION_VARIABLE m_hasWork;
	CONDITION_VARIABLE m_hasDone;
};


#define SMALL_FLOAT 0.001f

// #define LTR_DEBUG 1
#ifdef LTR_DEBUG
#define DBG( x ) x
#else
#define DBG( x )
#endif


#define LTR_WT_PREXFORM 0 // transform positions/normals for instances
#define LTR_WT_COLINFO  1 // generate collision info (* mesh)
#define LTR_WT_SAMPLES  2 // gather sampled information to compact position/normal/index arrays, allocate color arrays (* mesh)
#define LTR_WT_LMRENDER 3 // generate colors for samples by rendering lights to lightmaps (* mesh * light)
#define LTR_WT_RDGENLNK 4 // calculate sample links
#define LTR_WT_RDBOUNCE 5 // calculate a bounce between samples
#define LTR_WT_RDCOMMIT 6 // commit radiosity calculations
#define LTR_WT_AORENDER 7 // generate ambient occlusion for samples (* mesh)
#define LTR_WT_FINALIZE 8 // push sample data back to lightmaps (* mesh)


FORCEINLINE float randf(){ return (float) rand() / (float) RAND_MAX; }

double ltr_gettime();


template< class T > T TMIN( const T& a, const T& b ){ return a < b ? a : b; }
template< class T > T TMAX( const T& a, const T& b ){ return a > b ? a : b; }
template< class T > void TMEMSET( T* a, size_t c, const T& v )
{
	for( size_t i = 0; i < c; ++i )
		a[ i ] = v;
}
template< class T > void TMEMCOPY( T* a, const T* src, size_t c ){ memcpy( a, src, c * sizeof(T) ); }
template< class T, class S > T TLERP( const T& a, const T& b, const S& s ){ return a * ( S(1) - s ) + b * s; }


struct Vec2
{
	float x, y;
	
	static Vec2 Create( float x ){ Vec2 v = { x, x }; return v; }
	static Vec2 Create( float x, float y ){ Vec2 v = { x, y }; return v; }
	
	FORCEINLINE Vec2 operator + () const { return *this; }
	FORCEINLINE Vec2 operator - () const { Vec2 v = { -x, -y }; return v; }
	
	FORCEINLINE Vec2 operator + ( const Vec2& o ) const { Vec2 v = { x + o.x, y + o.y }; return v; }
	FORCEINLINE Vec2 operator - ( const Vec2& o ) const { Vec2 v = { x - o.x, y - o.y }; return v; }
	FORCEINLINE Vec2 operator * ( const Vec2& o ) const { Vec2 v = { x * o.x, y * o.y }; return v; }
	FORCEINLINE Vec2 operator / ( const Vec2& o ) const { Vec2 v = { x / o.x, y / o.y }; return v; }
	
	FORCEINLINE Vec2 operator + ( float f ) const { Vec2 v = { x + f, y + f }; return v; }
	FORCEINLINE Vec2 operator - ( float f ) const { Vec2 v = { x - f, y - f }; return v; }
	FORCEINLINE Vec2 operator * ( float f ) const { Vec2 v = { x * f, y * f }; return v; }
	FORCEINLINE Vec2 operator / ( float f ) const { Vec2 v = { x / f, y / f }; return v; }
	
	FORCEINLINE Vec2& operator += ( const Vec2& o ){ x += o.x; y += o.y; return *this; }
	FORCEINLINE Vec2& operator -= ( const Vec2& o ){ x -= o.x; y -= o.y; return *this; }
	FORCEINLINE Vec2& operator *= ( const Vec2& o ){ x *= o.x; y *= o.y; return *this; }
	FORCEINLINE Vec2& operator /= ( const Vec2& o ){ x /= o.x; y /= o.y; return *this; }
	
	FORCEINLINE Vec2& operator += ( float f ){ x += f; y += f; return *this; }
	FORCEINLINE Vec2& operator -= ( float f ){ x -= f; y -= f; return *this; }
	FORCEINLINE Vec2& operator *= ( float f ){ x *= f; y *= f; return *this; }
	FORCEINLINE Vec2& operator /= ( float f ){ x /= f; y /= f; return *this; }
	
	FORCEINLINE bool operator == ( const Vec2& o ) const { return x == o.x && y == o.y; }
	FORCEINLINE bool operator != ( const Vec2& o ) const { return x != o.x || y != o.y; }
	
	FORCEINLINE Vec2 Perp() const { Vec2 v = { y, -x }; return v; }
	FORCEINLINE float LengthSq() const { return x * x + y * y; }
	FORCEINLINE float Length() const { return sqrtf( LengthSq() ); }
	FORCEINLINE Vec2 Normalized() const
	{
		float lensq = LengthSq();
		if( lensq == 0 )
		{
			Vec2 v = { 0, 0 };
			return v;
		}
		float invlen = 1.0f / sqrtf( lensq );
		Vec2 v = { x * invlen, y * invlen };
		return v;
	}
	
	void Dump( FILE* f ) const
	{
		fprintf( f, "Vec2 ( %.2f %.2f )\n", x, y );
	}
};

FORCEINLINE float Vec2Dot( const Vec2& v1, const Vec2& v2 ){ return v1.x * v2.x + v1.y * v2.y; }
FORCEINLINE float Vec2Cross( const Vec2& v1, const Vec2& v2 )
{
	return ( v1.x * v2.y ) - ( v1.y * v2.x );
}

struct Vec3
{
	float x, y, z;
	
	static FORCEINLINE Vec3 Create( float x ){ Vec3 v = { x, x, x }; return v; }
	static FORCEINLINE Vec3 Create( float x, float y, float z ){ Vec3 v = { x, y, z }; return v; }
	static FORCEINLINE Vec3 CreateFromPtr( const float* x ){ Vec3 v = { x[0], x[1], x[2] }; return v; }
	static FORCEINLINE Vec3 CreateRandomVector( float maxdist )
	{
		float a = randf() * (float)M_PI * 2;
		float b = randf() * (float)M_PI;
		float d = randf() * maxdist;
		float ac = cos( a ), as = sin( a );
		float bc = cos( b ), bs = sin( b );
		Vec3 v = { ac * bs * d, as * bs * d, bc * d };
		return v;
	}
	static FORCEINLINE Vec3 CreateRandomVectorDirDvg( const Vec3& dir, float dvg );
	static FORCEINLINE Vec3 CreateSpiralDirVector( const Vec3& dir, float randoff, int i, int sample_count );
	
	FORCEINLINE Vec3 operator + () const { return *this; }
	FORCEINLINE Vec3 operator - () const { Vec3 v = { -x, -y, -z }; return v; }
	
	FORCEINLINE Vec3 operator + ( const Vec3& o ) const { Vec3 v = { x + o.x, y + o.y, z + o.z }; return v; }
	FORCEINLINE Vec3 operator - ( const Vec3& o ) const { Vec3 v = { x - o.x, y - o.y, z - o.z }; return v; }
	FORCEINLINE Vec3 operator * ( const Vec3& o ) const { Vec3 v = { x * o.x, y * o.y, z * o.z }; return v; }
	FORCEINLINE Vec3 operator / ( const Vec3& o ) const { Vec3 v = { x / o.x, y / o.y, z / o.z }; return v; }
	
	FORCEINLINE Vec3 operator + ( float f ) const { Vec3 v = { x + f, y + f, z + f }; return v; }
	FORCEINLINE Vec3 operator - ( float f ) const { Vec3 v = { x - f, y - f, z - f }; return v; }
	FORCEINLINE Vec3 operator * ( float f ) const { Vec3 v = { x * f, y * f, z * f }; return v; }
	FORCEINLINE Vec3 operator / ( float f ) const { Vec3 v = { x / f, y / f, z / f }; return v; }
	
	FORCEINLINE Vec3& operator += ( const Vec3& o ){ x += o.x; y += o.y; z += o.z; return *this; }
	FORCEINLINE Vec3& operator -= ( const Vec3& o ){ x -= o.x; y -= o.y; z -= o.z; return *this; }
	FORCEINLINE Vec3& operator *= ( const Vec3& o ){ x *= o.x; y *= o.y; z *= o.z; return *this; }
	FORCEINLINE Vec3& operator /= ( const Vec3& o ){ x /= o.x; y /= o.y; z /= o.z; return *this; }
	
	FORCEINLINE Vec3& operator += ( float f ){ x += f; y += f; z += f; return *this; }
	FORCEINLINE Vec3& operator -= ( float f ){ x -= f; y -= f; z -= f; return *this; }
	FORCEINLINE Vec3& operator *= ( float f ){ x *= f; y *= f; z *= f; return *this; }
	FORCEINLINE Vec3& operator /= ( float f ){ x /= f; y /= f; z /= f; return *this; }
	
	FORCEINLINE bool operator == ( const Vec3& o ) const { return x == o.x && y == o.y && z == o.z; }
	FORCEINLINE bool operator != ( const Vec3& o ) const { return x != o.x || y != o.y || z != o.z; }
	
	FORCEINLINE bool IsZero() const { return x == 0 && y == 0 && z == 0; }
	FORCEINLINE bool NearZero() const { return fabs(x) < SMALL_FLOAT && fabs(y) < SMALL_FLOAT && fabs(z) < SMALL_FLOAT; }
	FORCEINLINE float LengthSq() const { return x * x + y * y + z * z; }
	FORCEINLINE float Length() const { return sqrtf( LengthSq() ); }
	FORCEINLINE Vec3 Normalized() const
	{
		float lensq = LengthSq();
		if( lensq == 0 )
		{
			Vec3 v = { 0, 0, 0 };
			return v;
		}
		float invlen = 1.0f / sqrtf( lensq );
		Vec3 v = { x * invlen, y * invlen, z * invlen };
		return v;
	}
	
	void Dump( FILE* f ) const
	{
		fprintf( f, "Vec3 ( %.2f %.2f %.2f )\n", x, y, z );
	}
};

FORCEINLINE Vec3 operator + ( float f, const Vec3& v ){ Vec3 out = { f + v.x, f + v.y, f + v.z }; return out; }
FORCEINLINE Vec3 operator - ( float f, const Vec3& v ){ Vec3 out = { f - v.x, f - v.y, f - v.z }; return out; }
FORCEINLINE Vec3 operator * ( float f, const Vec3& v ){ Vec3 out = { f * v.x, f * v.y, f * v.z }; return out; }
FORCEINLINE Vec3 operator / ( float f, const Vec3& v ){ Vec3 out = { f / v.x, f / v.y, f / v.z }; return out; }

FORCEINLINE float Vec3Dot( const Vec3& v1, const Vec3& v2 ){ return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; }
FORCEINLINE Vec3 Vec3Cross( const Vec3& v1, const Vec3& v2 )
{
	Vec3 out =
	{
		v1.y * v2.z - v1.z * v2.y,
		v1.z * v2.x - v1.x * v2.z,
		v1.x * v2.y - v1.y * v2.x,
	};
	return out;
}

Vec3 Vec3::CreateRandomVectorDirDvg( const Vec3& dir, float dvg )
{
	float a = randf() * (float)M_PI * 2;
	float b = randf() * (float)M_PI * dvg;
	float ac = cos( a ), as = sin( a );
	float bc = cos( b ), bs = sin( b );
	Vec3 diffvec = { dir.y, -dir.z, dir.x };
	Vec3 up = Vec3Cross( dir, diffvec ).Normalized();
	Vec3 rt = Vec3Cross( dir, up );
	return ac * bs * rt + as * bs * up + bc * dir;
}

#define DEG2RAD( x ) ((x)/180.0f*(float)M_PI)
Vec3 Vec3::CreateSpiralDirVector( const Vec3& dir, float randoff, int i, int sample_count )
{
	float q = ( i + 0.5f ) / sample_count;
	float cos_side = sqrt( q );
	float sin_side = sin( acos( cos_side ) );
	float angle = ( i + randoff ) * DEG2RAD( 137.508f );
	float cos_around = cos( angle );
	float sin_around = sin( angle );
	
	Vec3 diffvec = { dir.y, -dir.z, dir.x };
	Vec3 up = Vec3Cross( dir, diffvec ).Normalized();
	Vec3 rt = Vec3Cross( dir, up );
	
	return cos_around * sin_side * rt + sin_around * sin_side * up + cos_side * dir;
}


struct Vec4
{
	float x, y, z, w;
	
	FORCEINLINE Vec4 operator + ( const Vec4& o ) const { Vec4 v = { x + o.x, y + o.y, z + o.z, w + o.w }; return v; }
	FORCEINLINE Vec4 operator - ( const Vec4& o ) const { Vec4 v = { x - o.x, y - o.y, z - o.z, w - o.w }; return v; }

	FORCEINLINE Vec4 operator * ( float f ) const { Vec4 v = { x * f, y * f, z * f, w * f }; return v; }
	
	Vec3 ToVec3() const { return Vec3::Create( x, y, z ); }
};
static FORCEINLINE Vec4 V4( float x ){ Vec4 o = { x, x, x, x }; return o; }
static FORCEINLINE Vec4 V4( float x, float y, float z, float w ){ Vec4 o = { x, y, z, w }; return o; }
static FORCEINLINE Vec4 V4( const Vec3& v, float w ){ Vec4 o = { v.x, v.y, v.z, w }; return o; }


struct Mat4
{
	union
	{
		float a[16];
		float m[4][4];
	};
	
	void SetIdentity()
	{
		for( int i = 0; i < 4; ++i )
			for( int j = 0; j < 4; ++j )
				m[i][j] = i == j;
	}
	static Mat4 CreateIdentity()
	{
		Mat4 m;
		m.SetIdentity();
		return m;
	}
	
	FORCEINLINE Vec3 Transform( const Vec3& v, float w ) const
	{
		Vec3 out =
		{
			v.x * m[0][0] + v.y * m[1][0] + v.z * m[2][0] + m[3][0] * w,
			v.x * m[0][1] + v.y * m[1][1] + v.z * m[2][1] + m[3][1] * w,
			v.x * m[0][2] + v.y * m[1][2] + v.z * m[2][2] + m[3][2] * w,
		};
		return out;
	}
	FORCEINLINE Vec3 TransformPos( const Vec3& pos ) const { return Transform( pos, 1.0f ); }
	FORCEINLINE Vec3 TransformNormal( const Vec3& nrm ) const { return Transform( nrm, 0.0f ); }
	
	bool InvertTo( Mat4& out );
	FORCEINLINE void Transpose()
	{
		std::swap( m[1][0], m[0][1] );
		std::swap( m[2][0], m[0][2] );
		std::swap( m[3][0], m[0][3] );
		std::swap( m[2][1], m[1][2] );
		std::swap( m[3][1], m[1][3] );
		std::swap( m[3][2], m[2][3] );
	}
	void GenNormalMatrix( Mat4& out ) const
	{
		out.m[0][0] = m[0][0]; out.m[0][1] = m[0][1]; out.m[0][2] = m[0][2];
		out.m[1][0] = m[1][0]; out.m[1][1] = m[1][1]; out.m[1][2] = m[1][2];
		out.m[2][0] = m[2][0]; out.m[2][1] = m[2][1]; out.m[2][2] = m[2][2];
		out.m[0][3] = out.m[1][3] = out.m[2][3] = 0;
		out.m[3][0] = out.m[3][1] = out.m[3][2] = 0;
		out.m[3][3] = 1;
		
		out.InvertTo( out );
		out.Transpose();
	}
	
	void Dump( FILE* f ) const
	{
		fprintf( f, "Mat4 (\n" );
		fprintf( f, "  %.2f %.2f %.2f %.2f\n", m[0][0], m[0][1], m[0][2], m[0][3] );
		fprintf( f, "  %.2f %.2f %.2f %.2f\n", m[1][0], m[1][1], m[1][2], m[1][3] );
		fprintf( f, "  %.2f %.2f %.2f %.2f\n", m[2][0], m[2][1], m[2][2], m[2][3] );
		fprintf( f, "  %.2f %.2f %.2f %.2f\n", m[3][0], m[3][1], m[3][2], m[3][3] );
		fprintf( f, ")\n" );
	}
};


typedef std::vector< u32 > U32Vector;
typedef std::vector< float > FloatVector;
typedef std::vector< Vec2 > Vec2Vector;
typedef std::vector< Vec3 > Vec3Vector;
typedef std::vector< Vec4 > Vec4Vector;
typedef std::vector< Mat4 > Mat4Vector;
typedef std::vector< ltr_WorkOutput > WorkOutputVector;


float TriangleArea( const Vec3& P1, const Vec3& P2, const Vec3& P3 );
float CalculateSampleArea( const Vec2& tex1, const Vec2& tex2, const Vec2& tex3, const Vec3& pos1, const Vec3& pos2, const Vec3& pos3 );

void TransformPositions( Vec3* out, Vec3* arr, size_t count, const Mat4& matrix );
void TransformNormals( Vec3* out, Vec3* arr, size_t count, const Mat4& matrix );
void RasterizeTriangle2D( Vec3* image, i32 width, i32 height, const Vec2& p1, const Vec2& p2, const Vec2& p3, const Vec3& v1, const Vec3& v2, const Vec3& v3 );
void RasterizeTriangle2D_x2_ex( Vec3* img1, Vec3* img2, Vec4* img3, i32 width, i32 height, float margin,
	const Vec2& p1, const Vec2& p2, const Vec2& p3,
	const Vec3& va1, const Vec3& va2, const Vec3& va3,
	const Vec3& vb1, const Vec3& vb2, const Vec3& vb3,
	const Vec4& vc1, const Vec4& vc2, const Vec4& vc3 );


void Generate_Gaussian_Kernel( float* out, int ext, float radius );
void Convolve_Transpose( float* src, float* dst, u32 width, u32 height, int blur_ext, float* kernel, float* tmp );


// BSP tree
// - best split plane is chosen by triangle normals and general direction of vertex positions (longest projection)
// - triangles are split to fit in the node

struct BSPTriangle
{
	Vec3 P1, P2, P3;
	
	bool CheckIsUseful() const
	{
		Vec3 e1 = P2 - P1, e2 = P3 - P1;
		return !Vec3Cross( e1, e2 ).NearZero();
	}
};
typedef std::vector< BSPTriangle > BSPTriVector;

struct BSPNode
{
	BSPNode() : front_node( NULL ), back_node( NULL ){}
	~BSPNode()
	{
		if( front_node ) delete front_node;
		if( back_node ) delete back_node;
	}
	
	void Split( int depth );
	void AddTriangleSplit( BSPTriangle* tri );
	float IntersectRay( const Vec3& from, const Vec3& to, Vec3* outnormal );
	bool PickSplitPlane();
	
	void Dump( FILE* f, int lev = 0, const char* pfx = "" )
	{
		for( int i = 0; i < lev; ++i )
			fputc( ' ', f );
		fprintf( f, "%sNODE [%g;%g;%g](%g), tris=%d\n", pfx, N.x, N.y, N.z, D, (int) triangles.size() );
		fprintf( f, "{\n" );
		for( size_t i = 0; i < triangles.size(); ++i )
		{
			fprintf( f, "  " ); triangles[i].P1.Dump( stdout );
			fprintf( f, "  " ); triangles[i].P2.Dump( stdout );
			fprintf( f, "  " ); triangles[i].P3.Dump( stdout );
		}
		if( front_node )
			front_node->Dump( f, lev + 1, "F " );
		if( back_node )
			back_node->Dump( f, lev + 1, "B " );
	}
	
	Vec3 N;
	float D;
	
	BSPNode *front_node, *back_node;
	
	BSPTriVector triangles;
};

struct BSPTree
{
	BSPTree() : root( new BSPNode() ){}
	~BSPTree(){ delete root; }
	
	FORCEINLINE void SetTriangles( BSPTriangle* tris, size_t count )
	{
		root->triangles.resize( count );
		TMEMCOPY( &root->triangles[0], tris, count );
		root->Split( 0 );
	}
	FORCEINLINE float IntersectRay( const Vec3& from, const Vec3& to, Vec3* outnormal = NULL ){ return root->IntersectRay( from, to, outnormal ); }
	
	BSPNode* root;
};


