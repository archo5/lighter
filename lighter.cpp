

#include "lighter_int.h"


struct ltr_MeshPart
{
	u32 m_vertexCount;
	u32 m_vertexOffset;
	u32 m_indexCount;
	u32 m_indexOffset;
};
typedef std::vector< ltr_MeshPart > MeshPartVector;

struct ltr_Mesh
{
	ltr_Mesh( ltr_Scene* s ) : m_scene( s ){}
	
	ltr_Scene* m_scene;
	std::string m_ident;
	Vec3Vector m_vpos;
	Vec3Vector m_vnrm;
	Vec2Vector m_vtex1;
	Vec2Vector m_vtex2;
	U32Vector m_indices;
	MeshPartVector m_parts;
};
typedef std::vector< ltr_Mesh* > MeshPtrVector;

struct ltr_MeshInstance
{
	// input
	ltr_Mesh* mesh;
	Mat4 matrix;
	u32 lm_width;
	u32 lm_height;
	
	// tmp
	Vec3Vector m_vpos;
	
	// output
	Vec3Vector m_samples_pos;
	Vec3Vector m_samples_nrm;
	U32Vector m_samples_loc;
	Vec3Vector m_lightmap;
};
typedef std::vector< ltr_MeshInstance > MeshInstVector;

struct ltr_Scene
{
	ltr_Scene() : m_workType( LTR_WT_COLINFO ), m_workPart( 0 ){}
	
	void DoWork();
	LTRCODE Advance();
	
	MeshPtrVector m_meshes;
	MeshInstVector m_meshInstances;
	LightVector m_lights;
	
	WorkOutputVector m_workOutput;
	
	u32 m_workType;
	u32 m_workPart;
};


ltr_Scene* ltr_CreateScene()
{
	return new ltr_Scene;
}

void ltr_DestroyScene( ltr_Scene* scene )
{
	delete scene;
}

void ltr_Scene::DoWork()
{
	switch( m_workType )
	{
	case LTR_WT_COLINFO:
		break;
	case LTR_WT_SAMPLES:
		break;
	case LTR_WT_LMRENDER:
		break;
	case LTR_WT_FINALIZE:
		break;
	}
}

LTRCODE ltr_Scene::Advance()
{
	if( m_workType == LTR_WT_FINALIZE && m_workPart == m_meshInstances.size() )
		return LTRC_ATEND;
	
	u32 count = m_meshInstances.size();
	if( m_workType == LTR_WT_LMRENDER )
		count *= m_lights.size();
	
	m_workPart++;
	if( m_workPart == count )
	{
		m_workType++;
		m_workPart = 0;
	}
	
	return LTRC_SUCCESS;
}

LTRCODE ltr_DoWork( ltr_Scene* scene, ltr_WorkInfo* info )
{
	scene->DoWork();
	
	info->type = scene->m_workType;
	info->part = scene->m_workPart;
	switch( scene->m_workType )
	{
	case LTR_WT_COLINFO:
		info->stage = "generating data structures";
		info->item_count = scene->m_meshInstances.size();
		break;
	case LTR_WT_SAMPLES:
		info->stage = "generating samples";
		info->item_count = scene->m_meshInstances.size();
		break;
	case LTR_WT_LMRENDER:
		info->stage = "rendering lightmaps";
		info->item_count = scene->m_meshInstances.size() * scene->m_lights.size();
		break;
	case LTR_WT_FINALIZE:
		info->stage = "exporting lightmaps";
		info->item_count = scene->m_meshInstances.size();
		break;
	default:
		info->stage = "unknown";
		info->item_count = 0;
		break;
	}
	
	return scene->Advance();
}


ltr_Mesh* ltr_CreateMesh( ltr_Scene* scene, const char* ident, size_t ident_size )
{
	ltr_Mesh* mesh = new ltr_Mesh( scene );
	if( ident )
		mesh->m_ident.assign( ident, ident_size );
	scene->m_meshes.push_back( mesh );
	return mesh;
}

void ltr_MeshAddPart( ltr_Mesh* mesh, ltr_MeshPartInfo* mpinfo )
{
	ltr_MeshPart mp = { mpinfo->vertex_count, (u32) mesh->m_vpos.size(), mpinfo->index_count, (u32) mesh->m_indices.size() };
	
	size_t newsize = mesh->m_vpos.size() + mpinfo->vertex_count;
	mesh->m_vpos.resize( newsize );
	mesh->m_vnrm.resize( newsize );
	mesh->m_vtex1.resize( newsize );
	mesh->m_vtex2.resize( newsize );
	mesh->m_indices.resize( mesh->m_indices.size() + mpinfo->index_count );
	
	Vec3* vpos_ptr = &mesh->m_vpos[ mp.m_vertexOffset ];
	Vec3* vnrm_ptr = &mesh->m_vnrm[ mp.m_vertexOffset ];
	Vec2* vtex1_ptr = &mesh->m_vtex1[ mp.m_vertexOffset ];
	Vec2* vtex2_ptr = &mesh->m_vtex2[ mp.m_vertexOffset ];
	u32* iptr = &mesh->m_indices[ mp.m_indexOffset ];
	
	for( u32 i = 0; i < mp.m_vertexCount; ++i )
	{
		vpos_ptr[ i ] = *(Vec3*)((char*) mpinfo->positions_f3 + mpinfo->stride_positions * i );
		vnrm_ptr[ i ] = *(Vec3*)((char*) mpinfo->normals_f3 + mpinfo->stride_normals * i );
		vtex1_ptr[ i ] = *(Vec2*)((char*) mpinfo->texcoords1_f2 + mpinfo->stride_texcoords1 * i );
		vtex2_ptr[ i ] = *(Vec2*)((char*) mpinfo->texcoords2_f2 + mpinfo->stride_texcoords2 * i );
	}
	memcpy( iptr, mpinfo->indices, sizeof(*iptr) * mpinfo->index_count );
	
	mesh->m_parts.push_back( mp );
}

LTRBOOL ltr_MeshAddInstance( ltr_Mesh* mesh, ltr_MeshInstanceInfo* mii )
{
	ltr_MeshInstance mi;
	mi.mesh = mesh;
	memcpy( mi.matrix.a, mii->matrix, sizeof(Mat4) );
	mesh->m_scene->m_meshInstances.push_back( mi );
	return 1;
}

void ltr_LightAdd( ltr_Scene* scene, ltr_LightInfo* li )
{
	scene->m_lights.push_back( *li );
}


void ltr_GetWorkOutputInfo( ltr_Scene* scene, ltr_WorkOutputInfo* woutinfo )
{
	woutinfo->lightmap_count = (i32) scene->m_workOutput.size();
}

LTRBOOL ltr_GetWorkOutput( ltr_Scene* scene, i32 which, ltr_WorkOutput* wout )
{
	if( which < 0 || which >= (i32) scene->m_workOutput.size() )
		return 0;
	*wout = scene->m_workOutput[ which ];
	return 1;
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
	float vcross = Vec2CrossProduct( p1p2, p1p3 );
	if( vcross == 0 )
		return;
	
	Vec3 va1va2 = va2 - va1, va1va3 = va3 - va1;
	
	for( i32 x = minX; x <= maxX; x++ )
	{
		for( i32 y = minY; y <= maxY; y++ )
		{
			Vec2 q = { x - p1.x, y - p1.y };
			
			float s = Vec2CrossProduct( q, p1p3 ) / vcross;
			float t = Vec2CrossProduct( p1p2, q ) / vcross;
			
			if( ( s >= 0 ) && ( t >= 0 ) && ( s + t <= 1 ) )
			{
				image[ x + width * y ] = va1 + va1va2 * s + va1va3 * t;
			}
		}
	}
}


