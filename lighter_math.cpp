

#include "lighter_int.hpp"


bool Mat4::InvertTo( Mat4& out )
{
	float inv[16], det;
	int i;
	
	inv[0] = a[5]  * a[10] * a[15] -
			 a[5]  * a[11] * a[14] -
			 a[9]  * a[6]  * a[15] +
			 a[9]  * a[7]  * a[14] +
			 a[13] * a[6]  * a[11] -
			 a[13] * a[7]  * a[10];
	
	inv[4] = -a[4]  * a[10] * a[15] +
			  a[4]  * a[11] * a[14] +
			  a[8]  * a[6]  * a[15] -
			  a[8]  * a[7]  * a[14] -
			  a[12] * a[6]  * a[11] +
			  a[12] * a[7]  * a[10];
	
	inv[8] = a[4]  * a[9] * a[15] -
			 a[4]  * a[11] * a[13] -
			 a[8]  * a[5] * a[15] +
			 a[8]  * a[7] * a[13] +
			 a[12] * a[5] * a[11] -
			 a[12] * a[7] * a[9];
	
	inv[12] = -a[4]  * a[9] * a[14] +
			   a[4]  * a[10] * a[13] +
			   a[8]  * a[5] * a[14] -
			   a[8]  * a[6] * a[13] -
			   a[12] * a[5] * a[10] +
			   a[12] * a[6] * a[9];
	
	inv[1] = -a[1]  * a[10] * a[15] +
			  a[1]  * a[11] * a[14] +
			  a[9]  * a[2] * a[15] -
			  a[9]  * a[3] * a[14] -
			  a[13] * a[2] * a[11] +
			  a[13] * a[3] * a[10];
	
	inv[5] = a[0]  * a[10] * a[15] -
			 a[0]  * a[11] * a[14] -
			 a[8]  * a[2] * a[15] +
			 a[8]  * a[3] * a[14] +
			 a[12] * a[2] * a[11] -
			 a[12] * a[3] * a[10];
	
	inv[9] = -a[0]  * a[9] * a[15] +
			  a[0]  * a[11] * a[13] +
			  a[8]  * a[1] * a[15] -
			  a[8]  * a[3] * a[13] -
			  a[12] * a[1] * a[11] +
			  a[12] * a[3] * a[9];
	
	inv[13] = a[0]  * a[9] * a[14] -
			  a[0]  * a[10] * a[13] -
			  a[8]  * a[1] * a[14] +
			  a[8]  * a[2] * a[13] +
			  a[12] * a[1] * a[10] -
			  a[12] * a[2] * a[9];
	
	inv[2] = a[1]  * a[6] * a[15] -
			 a[1]  * a[7] * a[14] -
			 a[5]  * a[2] * a[15] +
			 a[5]  * a[3] * a[14] +
			 a[13] * a[2] * a[7] -
			 a[13] * a[3] * a[6];
	
	inv[6] = -a[0]  * a[6] * a[15] +
			  a[0]  * a[7] * a[14] +
			  a[4]  * a[2] * a[15] -
			  a[4]  * a[3] * a[14] -
			  a[12] * a[2] * a[7] +
			  a[12] * a[3] * a[6];
	
	inv[10] = a[0]  * a[5] * a[15] -
			  a[0]  * a[7] * a[13] -
			  a[4]  * a[1] * a[15] +
			  a[4]  * a[3] * a[13] +
			  a[12] * a[1] * a[7] -
			  a[12] * a[3] * a[5];
	
	inv[14] = -a[0]  * a[5] * a[14] +
			   a[0]  * a[6] * a[13] +
			   a[4]  * a[1] * a[14] -
			   a[4]  * a[2] * a[13] -
			   a[12] * a[1] * a[6] +
			   a[12] * a[2] * a[5];
	
	inv[3] = -a[1] * a[6] * a[11] +
			  a[1] * a[7] * a[10] +
			  a[5] * a[2] * a[11] -
			  a[5] * a[3] * a[10] -
			  a[9] * a[2] * a[7] +
			  a[9] * a[3] * a[6];
	
	inv[7] = a[0] * a[6] * a[11] -
			 a[0] * a[7] * a[10] -
			 a[4] * a[2] * a[11] +
			 a[4] * a[3] * a[10] +
			 a[8] * a[2] * a[7] -
			 a[8] * a[3] * a[6];
	
	inv[11] = -a[0] * a[5] * a[11] +
			   a[0] * a[7] * a[9] +
			   a[4] * a[1] * a[11] -
			   a[4] * a[3] * a[9] -
			   a[8] * a[1] * a[7] +
			   a[8] * a[3] * a[5];
	
	inv[15] = a[0] * a[5] * a[10] -
			  a[0] * a[6] * a[9] -
			  a[4] * a[1] * a[10] +
			  a[4] * a[2] * a[9] +
			  a[8] * a[1] * a[6] -
			  a[8] * a[2] * a[5];
	
	det = a[0] * inv[0] + a[1] * inv[4] + a[2] * inv[8] + a[3] * inv[12];
	
	if( det == 0 )
		return false;
	
	det = 1.0f / det;
	
	for( i = 0; i < 16; ++i )
		out.a[ i ] = inv[ i ] * det;
	
	return true;
}


static FORCEINLINE float TriangleArea( float a, float b, float c )
{
	float p = ( a + b + c ) * 0.5f;
	float presqrt = p * ( p - a ) * ( p - b ) * ( p - c );
	if( presqrt < 0 )
		return 0;
	return sqrtf( presqrt );
}

float TriangleArea( const Vec2& P1, const Vec2& P2, const Vec2& P3 )
{
	float a = ( P2 - P1 ).Length();
	float b = ( P3 - P2 ).Length();
	float c = ( P1 - P3 ).Length();
	return TriangleArea( a, b, c );
}

float TriangleArea( const Vec3& P1, const Vec3& P2, const Vec3& P3 )
{
	float a = ( P2 - P1 ).Length();
	float b = ( P3 - P2 ).Length();
	float c = ( P1 - P3 ).Length();
	return TriangleArea( a, b, c );
}

float CalculateSampleArea( const Vec2& tex1, const Vec2& tex2, const Vec2& tex3, const Vec3& pos1, const Vec3& pos2, const Vec3& pos3 )
{
	float lmarea = TriangleArea( tex1, tex2, tex3 );
	float coarea = TriangleArea( pos1, pos2, pos3 );
	
	return lmarea > 0 ? coarea / lmarea : 0;
}


//
// TRANSFORMATION
//
void TransformPositions( Vec3* out, Vec3* arr, size_t count, const Mat4& matrix )
{
	for( size_t i = 0; i < count; ++i )
		out[i] = matrix.TransformPos( arr[i] );
}
void TransformNormals( Vec3* out, Vec3* arr, size_t count, const Mat4& matrix )
{
	Mat4 nrmtx;
	matrix.GenNormalMatrix( nrmtx );
	for( size_t i = 0; i < count; ++i )
		out[i] = nrmtx.TransformNormal( arr[i] ).Normalized();
}


//
// RASTERIZATION
//
// - p1, p2, p3 - in image space
//
void RasterizeTriangle2D( Vec3* image, i32 width, i32 height, const Vec2& p1, const Vec2& p2, const Vec2& p3, const Vec3& va1, const Vec3& va2, const Vec3& va3 )
{
	i32 maxX = (i32) TMAX( p1.x, TMAX( p2.x, p3.x ) );
	i32 minX = (i32) TMIN( p1.x, TMIN( p2.x, p3.x ) );
	i32 maxY = (i32) TMAX( p1.y, TMAX( p2.y, p3.y ) );
	i32 minY = (i32) TMIN( p1.y, TMIN( p2.y, p3.y ) );
	
	if( maxX < 0 || minX >= width || maxY < 0 || minY >= height )
		return;
	if( maxX < 0 ) maxX = 0;
	if( minX >= width ) minX = width - 1;
	if( maxY < 0 ) maxY = 0;
	if( minY >= height ) minY = height - 1;
	
	Vec2 p1p2 = p2 - p1;
	Vec2 p1p3 = p3 - p1;
	float vcross = Vec2Cross( p1p2, p1p3 );
	if( vcross == 0 )
		return;
	
	Vec3 va1va2 = va2 - va1, va1va3 = va3 - va1;
	
	for( i32 x = minX; x <= maxX; x++ )
	{
		for( i32 y = minY; y <= maxY; y++ )
		{
			Vec2 q = { x - p1.x, y - p1.y };
			
			float s = Vec2Cross( q, p1p3 ) / vcross;
			float t = Vec2Cross( p1p2, q ) / vcross;
			
			if( ( s >= 0 ) && ( t >= 0 ) && ( s + t <= 1 ) )
			{
				image[ x + width * y ] = va1 + va1va2 * s + va1va3 * t;
			}
		}
	}
}

void RasterizeTriangle2D_x2_ex( Vec3* img1, Vec3* img2, Vec4* img3, i32 width, i32 height, float margin,
	const Vec2& p1, const Vec2& p2, const Vec2& p3,
	const Vec3& va1, const Vec3& va2, const Vec3& va3,
	const Vec3& vb1, const Vec3& vb2, const Vec3& vb3,
	const Vec4& vc1, const Vec4& vc2, const Vec4& vc3 )
{
	i32 maxX = i32( TMAX( p1.x, TMAX( p2.x, p3.x ) ) + margin );
	i32 minX = i32( TMIN( p1.x, TMIN( p2.x, p3.x ) ) - margin );
	i32 maxY = i32( TMAX( p1.y, TMAX( p2.y, p3.y ) ) + margin );
	i32 minY = i32( TMIN( p1.y, TMIN( p2.y, p3.y ) ) - margin );
	
	if( maxX < 0 || minX >= width || maxY < 0 || minY >= height )
		return;
	if( minX < 0 ) minX = 0;
	if( maxX >= width ) maxX = width - 1;
	if( minY < 0 ) minY = 0;
	if( maxY >= height ) maxY = height - 1;
	
	Vec2 p1p2 = p2 - p1;
	Vec2 p1p3 = p3 - p1;
	float vcross = Vec2Cross( p1p2, p1p3 );
	if( vcross == 0 )
		return;
	
	Vec2 n1 = p1p2.Perp().Normalized();
	Vec2 n2 = ( p3 - p2 ).Perp().Normalized();
	Vec2 n3 = -p1p3.Perp().Normalized();
	float d1 = Vec2Dot( n1, p1 );
	float d2 = Vec2Dot( n2, p2 );
	float d3 = Vec2Dot( n3, p3 );
	float MG = margin;
	
	Vec3 va1va2 = va2 - va1, va1va3 = va3 - va1;
	Vec3 vb1vb2 = vb2 - vb1, vb1vb3 = vb3 - vb1;
	Vec4 vc1vc2 = vc2 - vc1, vc1vc3 = vc3 - vc1;
	
	for( i32 x = minX; x <= maxX; x++ )
	{
		for( i32 y = minY; y <= maxY; y++ )
		{
			Vec2 q = { x - p1.x, y - p1.y };
			
			float s = Vec2Cross( q, p1p3 ) / vcross;
			float t = Vec2Cross( p1p2, q ) / vcross;
			
			Vec2 p = { (float) x, (float) y };
			float pd1 = Vec2Dot( p, n1 ), pd2 = Vec2Dot( p, n2 ), pd3 = Vec2Dot( p, n3 );
			if( ( pd1 <= d1 + MG && pd2 <= d2 + MG && pd3 <= d3 + MG ) || ( pd1 + MG >= d1 && pd2 + MG >= d2 && pd3 + MG >= d3 ) )
			{
				img1[ x + width * y ] = va1 + va1va2 * s + va1va3 * t;
				img2[ x + width * y ] = vb1 + vb1vb2 * s + vb1vb3 * t;
				img3[ x + width * y ] = vc1 + vc1vc2 * s + vc1vc3 * t;
			}
		}
	}
}


float PlaneTriangleIntersect( const Vec3& N, float D, const Vec3& P1, const Vec3& P2, const Vec3& P3 )
{
	float proj1 = Vec3Dot( N, P1 ) - D;
	float proj2 = Vec3Dot( N, P2 ) - D;
	float proj3 = Vec3Dot( N, P3 ) - D;
	if( proj1 * proj2 > 0 && proj1 * proj3 > 0 )
		return proj1;
	return 0;
}

float IntersectLineSegmentTriangle( const Vec3& L1, const Vec3& L2, const Vec3& P1, const Vec3& P2, const Vec3& P3 )
{
	Vec3 u = P2 - P1;
	Vec3 v = P3 - P1;
	Vec3 n = Vec3Cross( u, v );
	if( n.NearZero() )
		return 2.0f;
	
	Vec3 dir = L2 - L1;
	Vec3 w0 = L1 - P1;
	float a = -Vec3Dot( n, w0 );
	float b = Vec3Dot( n, dir );
	if( fabs( b ) < SMALL_FLOAT )
		return 2.0f;
	
	float r = a / b;
	if( r < 0.0f || r > 1.0f )
	    return 2.0f;
	
	Vec3 I = L1 + r * dir;
	
	// is I inside T?
	float uu = Vec3Dot( u, u );
	float uv = Vec3Dot( u, v );
	float vv = Vec3Dot( v, v );
	Vec3 w = I - P1;
	float wu = Vec3Dot( w, u );
	float wv = Vec3Dot( w, v );
	float D = uv * uv - uu * vv;
	
	// get and test parametric coords
	float s = ( uv * wv - vv * wu ) / D;
	if( s < 0.0f || s > 1.0f )
		return 2.0f;
	float t = ( uv * wu - uu * wv ) / D;
	if( t < 0.0f || ( s + t ) > 1.0f )
		return 2.0f;
	
	return r;
}


void Generate_Gaussian_Kernel( float* out, int ext, float radius )
{
	float sum = 0.0f;
	float multconst = 1.0f / sqrtf( 2.0f * (float) M_PI * radius * radius );
	for( int i = -ext; i <= ext; ++i )
		sum += out[ i + ext ] = exp( -0.5f * pow( (float) i / radius, 2.0f ) ) * multconst;
	for( int i = -ext; i <= ext; ++i )
		out[ i + ext ] /= sum;
}

void Convolve_Transpose( float* src, float* dst, u32 width, u32 height, int conext, float* kernel, float* tmp )
{
	for( u32 y = 0; y < height; ++y )
	{
		// write row to temporary buffer, extend data beyond sides
		memcpy( tmp + conext * 3, src + y * width * 3, sizeof(float) * width * 3 );
		for( int i = 0; i < conext; ++i )
		{
			tmp[ i * 3 + 0 ] = tmp[ conext * 3 + 0 ];
			tmp[ i * 3 + 1 ] = tmp[ conext * 3 + 1 ];
			tmp[ i * 3 + 2 ] = tmp[ conext * 3 + 2 ];
			tmp[ ( conext + width + i ) * 3 + 0 ] = tmp[ ( conext + width - 1 ) * 3 + 0 ];
			tmp[ ( conext + width + i ) * 3 + 1 ] = tmp[ ( conext + width - 1 ) * 3 + 1 ];
			tmp[ ( conext + width + i ) * 3 + 2 ] = tmp[ ( conext + width - 1 ) * 3 + 2 ];
		}
		// convolve with the kernel
		for( u32 x = 0; x < width; ++x )
		{
			float* src_p = tmp + ( x + conext ) * 3;
			float* dst_p = dst + ( y + height * x ) * 3;
			float sum[3] = { 0.0f, 0.0f, 0.0f };
			for( int i = -conext; i <= conext; ++i )
			{
				float kq = kernel[ i + conext ];
				sum[0] += src_p[ i * 3 + 0 ] * kq;
				sum[1] += src_p[ i * 3 + 1 ] * kq;
				sum[2] += src_p[ i * 3 + 2 ] * kq;
			}
			LTR_VEC3_SET( dst_p, sum[0], sum[1], sum[2] );
		}
	}
}


#define BSP_NO_HIT 2.0f
#define BSP_MAX_NODE_COUNT 16
#define BSP_MAX_NODE_DEPTH 32
#define BSP_REEVAL_INTERVAL 8

void BSPNode::AddTriangle( BSPTriangle* tri, int depth )
{
	if( !tri->CheckIsUseful() )
		return;
	
	if( front_node ) // node already split
	{
		AddTriangleSplit( tri, depth + 1 );
	}
	else
	{
		triangles.push_back( *tri );
		if( triangles.size() > BSP_MAX_NODE_COUNT &&
			triangles.size() % BSP_REEVAL_INTERVAL == 0 &&
			depth < BSP_MAX_NODE_DEPTH )
		{
			if( PickSplitPlane() )
			{
				front_node = new BSPNode;
				back_node = new BSPNode;
				for( size_t i = 0; i < triangles.size(); ++i )
					AddTriangleSplit( &triangles[i], depth + 1 );
				BSPTriVector().swap( triangles );
			}
		}
	}
}

void BSPNode::AddTriangleSplit( BSPTriangle* tri, int depth )
{
	Vec3 P1 = tri->P1;
	Vec3 P2 = tri->P2;
	Vec3 P3 = tri->P3;
	
	// assume that plane is picked (valid N, D)
	float proj1 = Vec3Dot( N, P1 ) - D;
	float proj2 = Vec3Dot( N, P2 ) - D;
	float proj3 = Vec3Dot( N, P3 ) - D;
	if( proj1 * proj2 >= -SMALL_FLOAT && proj1 * proj3 >= -SMALL_FLOAT )
	{
		( proj1 + proj2 + proj3 > 0 ? front_node : back_node )->AddTriangle( tri, depth );
		return;
	}
	
	Vec3 S1 = P1, S2 = P2, S3 = P3;
	float td12 = proj1 - proj2; // opposite signs expected, both must be 0 for diff to be 0
	if( td12 )
		S1 = TLERP( P1, P2, fabs( proj1 / td12 ) );
	float td23 = proj2 - proj3; // opposite signs expected, both must be 0 for diff to be 0
	if( td23 )
		S2 = TLERP( P2, P3, fabs( proj2 / td23 ) );
	float td31 = proj3 - proj1; // opposite signs expected, both must be 0 for diff to be 0
	if( td31 )
		S3 = TLERP( P3, P1, fabs( proj3 / td31 ) );
	
	static int wat = 0;
	if( wat++ > 10000 )
		puts( "TODO FIX LEAK" );
	
	front_node->AddTriangle( tri, depth );
	back_node->AddTriangle( tri, depth );
}

float BSPNode::IntersectRay( const Vec3& from, const Vec3& to, Vec3* outnormal )
{
	if( front_node ) // node split
	{
		float d_from = Vec3Dot( from, N ) - D;
		float d_to = Vec3Dot( to, N ) - D;
		if( d_from < 0 )
		{
			float hit = back_node->IntersectRay( from, to, outnormal );
			if( hit < 1.0f )
				return hit;
			return d_to > 0 ? front_node->IntersectRay( from, to, outnormal ) : BSP_NO_HIT;
		}
		else
		{
			float hit = front_node->IntersectRay( from, to, outnormal );
			if( hit < 1.0f )
				return hit;
			return d_to < 0 ? back_node->IntersectRay( from, to, outnormal ) : BSP_NO_HIT;
		}
	}
	else
	{
		float closest_hit = BSP_NO_HIT;
		BSPTriangle* closest_tri = NULL;
		for( size_t i = 0; i < triangles.size(); ++i )
		{
			BSPTriangle& T = triangles[i];
			float hit = IntersectLineSegmentTriangle( from, to, T.P1, T.P2, T.P3 );
			if( hit < closest_hit )
			{
				closest_hit = hit;
				closest_tri = &T;
				// printf( "HIT %f from=%.2f %.2f %.2f to=%.2f %.2f %.2f p1=%.4f %.4f %.4f p2=%.4f %.4f %.4f p3=%.4f %.4f %.4f",
				// hit,from.x,from.y,from.z,to.x,to.y,to.z,T.P1.x,T.P1.y,T.P1.z,T.P2.x,T.P2.y,T.P2.z,T.P3.x,T.P3.y,T.P3.z);
			}
		}
		if( outnormal && closest_tri )
		{
			*outnormal = Vec3Cross( closest_tri->P3 - closest_tri->P1, closest_tri->P2 - closest_tri->P1 ).Normalized();
		}
		return closest_hit;
	}
}

bool BSPNode::PickSplitPlane()
{
	// split direction by positions / normals
	float mult = 1.0f / triangles.size();
	float mult3 = mult / 3.0f;
	Vec3 PN = {0,0,0}, NN = {0,0,0};
	for( size_t i = 0; i < triangles.size(); ++i )
	{
		BSPTriangle& T = triangles[i];
		Vec3 ND1 = ( T.P2 - T.P1 ) * mult3;
		Vec3 ND2 = ( T.P3 - T.P2 ) * mult3;
		Vec3 ND3 = ( T.P1 - T.P3 ) * mult3;
		Vec3 NNRM = Vec3Cross( ND1, ND3 ).Normalized();
		PN += Vec3Dot( PN, ND1 ) < 0 ? -ND1 : ND1;
		PN += Vec3Dot( PN, ND2 ) < 0 ? -ND2 : ND2;
		PN += Vec3Dot( PN, ND3 ) < 0 ? -ND3 : ND3;
		NN += Vec3Dot( NN, NNRM ) < 0 ? -NNRM : NNRM;
	}
	N = ( PN.Normalized() * 0.5f + NN.Normalized() * 0.5f ).Normalized();
	
	float dsum = 0;
	for( size_t i = 0; i < triangles.size(); ++i )
	{
		BSPTriangle& T = triangles[i];
		dsum += Vec3Dot( T.P1, N ) + Vec3Dot( T.P2, N ) + Vec3Dot( T.P3, N );
	}
	D = dsum * mult3;
	
	// evaluate viability of splitting
	int numA = 0;
	int numB = 0;
	float d1, d2, d3;
	for( size_t i = 0; i < triangles.size(); ++i )
	{
		BSPTriangle& T = triangles[i];
		d1 = Vec3Dot( T.P1, N ) - D;
		d2 = Vec3Dot( T.P2, N ) - D;
		d3 = Vec3Dot( T.P3, N ) - D;
		if( d1 * d2 > -SMALL_FLOAT && d2 * d3 > -SMALL_FLOAT )
		{
			if( d1 > 0 )
				numA++;
			else
				numB++;
		}
		else
		{
			numA += ( d1 > 0 ) + ( d2 > 0 ) + ( d3 > 0 );
			numB += ( d1 <=0 ) + ( d2 <=0 ) + ( d3 <=0 );
		}
	}
	if( numA != 0 && numB != 0 && numA + numB < triangles.size() * 1.7f )
	{
		float q = (float) numA / (float) numB;
		if( q < 1 ) q = 1 / q;
		return q < 3.0f;
	}
	return false;
}


