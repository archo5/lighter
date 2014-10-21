

#include "lighter_int.hpp"


void Mat4::InvertTo( Mat4& out )
{
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
		out[i] = nrmtx.TransformNormal( arr[i] );
}


//
// RASTERIZATION
//
// - p1, p2, p3 - in image space
//
void RasterizeTriangle2D( Vec3* image, i32 width, i32 height, const Vec2& p1, const Vec2& p2, const Vec2& p3, const Vec3& va1, const Vec3& va2, const Vec3& va3 )
{
	i32 maxX = TMAX( p1.x, TMAX( p2.x, p3.x ) );
	i32 minX = TMIN( p1.x, TMIN( p2.x, p3.x ) );
	i32 maxY = TMAX( p1.y, TMAX( p2.y, p3.y ) );
	i32 minY = TMIN( p1.y, TMIN( p2.y, p3.y ) );
	
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

void RasterizeTriangle2D_x2_ex( Vec3* img1, Vec3* img2, i32 width, i32 height, float margin,
	const Vec2& p1, const Vec2& p2, const Vec2& p3,
	const Vec3& va1, const Vec3& va2, const Vec3& va3,
	const Vec3& vb1, const Vec3& vb2, const Vec3& vb3 )
{
	i32 maxX = TMAX( p1.x, TMAX( p2.x, p3.x ) ) + margin;
	i32 minX = TMIN( p1.x, TMIN( p2.x, p3.x ) ) - margin;
	i32 maxY = TMAX( p1.y, TMAX( p2.y, p3.y ) ) + margin;
	i32 minY = TMIN( p1.y, TMIN( p2.y, p3.y ) ) - margin;
	
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
	
	float margin_ep12 = margin / p1p2.Length();
	float margin_ep13 = margin / p1p3.Length();
	float marginsum_1 = margin_ep12 + margin_ep13 + 1;
	
	Vec3 va1va2 = va2 - va1, va1va3 = va3 - va1;
	Vec3 vb1vb2 = vb2 - vb1, vb1vb3 = vb3 - vb1;
	
	for( i32 x = minX; x <= maxX; x++ )
	{
		for( i32 y = minY; y <= maxY; y++ )
		{
			Vec2 q = { x - p1.x, y - p1.y };
			
			float s = Vec2Cross( q, p1p3 ) / vcross;
			float t = Vec2Cross( p1p2, q ) / vcross;
			
			if( ( s >= -margin_ep13 ) && ( t >= -margin_ep12 ) && ( s + t <= marginsum_1 ) )
			{
				img1[ x + width * y ] = va1 + va1va2 * s + va1va3 * t;
				img2[ x + width * y ] = vb1 + vb1vb2 * s + vb1vb3 * t;
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
	Vec3 E0 = P2 - P1, E1 = P3 - P1;
	if( E0.IsZero() || E1.IsZero() )
		return 2.0f;
	Vec3 N = Vec3Cross( E0, E1 ).Normalized();
	float D = Vec3Dot( N, P1 ), D0 = Vec3Dot( N, L1 ), D1 = Vec3Dot( N, L2 );
	if( ( D0 - D ) * ( D1 - D ) >= 0 )
		return 2.0f;
	float dist = ( D - D0 ) / ( D1 - D0 );
	
	Vec3 isp = L1 * ( 1 - dist ) + L2 * dist;
	float s = ( Vec3Dot( isp, E0 ) - Vec3Dot( P1, E0 ) ) / E0.LengthSq();
	float t = ( Vec3Dot( isp, E1 ) - Vec3Dot( P1, E1 ) ) / E1.LengthSq();
	if( s <= 0 || t <= 0 || s + t >= 1 )
		return 2.0f;
	return dist;
}


#define BSP_NO_HIT 2.0f
#define BSP_MAX_NODE_COUNT 16
#define BSP_MAX_NODE_DEPTH 32

void BSPNode::AddTriangle( BSPTriangle* tri, int depth )
{
	if( front_node ) // node already split
	{
		float res = PlaneTriangleIntersect( N, D, tri->P1, tri->P2, tri->P3 );
		if( res == 0 )
		{
			AddTriangleSplit( tri, depth + 1 );
		}
		else
		{
			( res > 0 ? front_node : back_node )->AddTriangle( tri, depth + 1 );
		}
	}
	else
	{
		triangles.push_back( *tri );
		if( triangles.size() > BSP_MAX_NODE_COUNT && depth < BSP_MAX_NODE_DEPTH )
		{
			PickSplitPlane();
			for( size_t i = 0; i < triangles.size(); ++i )
				AddTriangleSplit( tri, depth + 1 );
			std::vector< BSPTriangle >().swap( triangles );
		}
	}
}

void BSPNode::AddTriangleSplit( BSPTriangle* tri, int depth )
{
	// TODO
}

float BSPNode::IntersectRay( const Vec3& from, const Vec3& to )
{
	if( front_node ) // node already split
	{
		float d_from = Vec3Dot( from, N ) - D;
		float d_to = Vec3Dot( to, N ) - D;
		if( d_from < 0 )
		{
			float hit = back_node->IntersectRay( from, to );
			if( hit < 1.0f )
				return hit;
			return d_to > 0 ? front_node->IntersectRay( from, to ) : BSP_NO_HIT;
		}
		else
		{
			float hit = front_node->IntersectRay( from, to );
			if( hit < 1.0f )
				return hit;
			return d_to < 0 ? back_node->IntersectRay( from, to ) : BSP_NO_HIT;
		}
	}
	else
	{
		float closest_hit = BSP_NO_HIT;
		for( size_t i = 0; i < triangles.size(); ++i )
		{
			BSPTriangle& T = triangles[i];
			float hit = IntersectLineSegmentTriangle( from, to, T.P1, T.P2, T.P3 );
			if( hit < closest_hit )
				closest_hit = hit;
		}
		return closest_hit;
	}
}

void BSPNode::PickSplitPlane()
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
		PN += Vec3Dot( PN, ND1 ) < 0 ? -ND1 : ND1;
		PN += Vec3Dot( PN, ND2 ) < 0 ? -ND2 : ND2;
		PN += Vec3Dot( PN, ND3 ) < 0 ? -ND3 : ND3;
		NN += Vec3Cross( ND1, -ND3 ).Normalized() * mult;
	}
	N = ( PN.Normalized() * 0.25f + NN.Normalized() * 0.75f ).Normalized();
	
	float dsum = 0;
	for( size_t i = 0; i < triangles.size(); ++i )
	{
		BSPTriangle& T = triangles[i];
		dsum += Vec3Dot( T.P1, N ) + Vec3Dot( T.P2, N ) + Vec3Dot( T.P3, N );
	}
	D = dsum * mult3;
}


