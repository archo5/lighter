

#include "lighter_int.hpp"


struct ltr_MeshPart
{
	u32 m_vertexCount;
	u32 m_vertexOffset;
	u32 m_indexCount;
	u32 m_indexOffset;
	int m_tristrip;
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
	Vec3Vector m_vnrm;
	Vec2Vector m_ltex;
	
	// output
	Vec3Vector m_samples_pos;
	Vec3Vector m_samples_nrm;
	U32Vector m_samples_loc;
	Vec3Vector m_lightmap;
};
typedef std::vector< ltr_MeshInstance* > MeshInstPtrVector;

struct ltr_Scene
{
	ltr_Scene() : m_workType( LTR_WT_PREXFORM ), m_workPart( 0 ){}
	~ltr_Scene()
	{
		for( size_t i = 0; i < m_meshInstances.size(); ++i )
			delete m_meshInstances[i];
		for( size_t i = 0; i < m_meshes.size(); ++i )
			delete m_meshes[i];
	}
	
	void DoWork();
	LTRCODE Advance();
	void RasterizeInstance( ltr_MeshInstance* mi, float margin );
	
	MeshPtrVector m_meshes;
	MeshInstPtrVector m_meshInstances;
	LightVector m_lights;
	
	Vec3Vector m_tmpRender1;
	Vec3Vector m_tmpRender2;
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

void ltr_Scene::RasterizeInstance( ltr_MeshInstance* mi, float margin )
{
	ltr_Mesh* mesh = mi->mesh;
	for( u32 part = 0; part < mesh->m_parts.size(); ++part )
	{
		const ltr_MeshPart& mp = mesh->m_parts[ part ];
		if( mp.m_tristrip )
		{
			for( u32 tri = 2; tri < mp.m_indexCount; ++tri )
			{
				u32 tridx1 = tri, tridx2 = tri + 1 + tri % 2, tridx3 = tri + 2 - tri % 2;
				RasterizeTriangle2D_x2_ex
				(
					&m_tmpRender1[0], &m_tmpRender2[0], mi->lm_width, mi->lm_height, margin,
					mi->m_ltex[ tridx1 ], mi->m_ltex[ tridx2 ], mi->m_ltex[ tridx3 ],
					mi->m_vpos[ tridx1 ], mi->m_vpos[ tridx2 ], mi->m_vpos[ tridx3 ],
					mi->m_vnrm[ tridx1 ], mi->m_vnrm[ tridx2 ], mi->m_vnrm[ tridx3 ]
				);
			}
		}
		else
		{
			for( u32 tri = 0; tri < mp.m_indexCount; tri += 3 )
			{
				u32 tridx1 = tri, tridx2 = tri + 1, tridx3 = tri + 2;
				RasterizeTriangle2D_x2_ex
				(
					&m_tmpRender1[0], &m_tmpRender2[0], mi->lm_width, mi->lm_height, margin,
					mi->m_ltex[ tridx1 ], mi->m_ltex[ tridx2 ], mi->m_ltex[ tridx3 ],
					mi->m_vpos[ tridx1 ], mi->m_vpos[ tridx2 ], mi->m_vpos[ tridx3 ],
					mi->m_vnrm[ tridx1 ], mi->m_vnrm[ tridx2 ], mi->m_vnrm[ tridx3 ]
				);
			}
		}
	}
}

void ltr_Scene::DoWork()
{
	switch( m_workType )
	{
	case LTR_WT_PREXFORM:
		{
			ltr_MeshInstance* mi = m_meshInstances[ m_workPart ];
			ltr_Mesh* mesh = mi->mesh;
			
			mi->m_vpos.resize( mesh->m_vpos.size() );
			mi->m_vnrm.resize( mesh->m_vnrm.size() );
			mi->m_ltex.resize( mesh->m_vtex2.size() );
			
			TransformPositions( &mi->m_vpos[0], &mesh->m_vpos[0], mesh->m_vpos.size(), mi->matrix );
			TransformNormals( &mi->m_vnrm[0], &mesh->m_vnrm[0], mesh->m_vnrm.size(), mi->matrix );
			
			Vec2 lsize = Vec2::Create( mi->lm_width, mi->lm_height );
			for( size_t i = 0; i < mi->m_ltex.size(); ++i )
				mi->m_ltex[ i ] = mesh->m_vtex2[ i ] * lsize;
		}
		break;
		
	case LTR_WT_COLINFO:
		break;
		
	case LTR_WT_SAMPLES:
		{
			ltr_MeshInstance* mi = m_meshInstances[ m_workPart ];
			
			if( m_tmpRender1.size() < mi->lm_width * mi->lm_height )
				m_tmpRender1.resize( mi->lm_width * mi->lm_height );
			if( m_tmpRender2.size() < mi->lm_width * mi->lm_height )
				m_tmpRender2.resize( mi->lm_width * mi->lm_height );
			
			TMEMSET( &m_tmpRender1[0], m_tmpRender1.size(), Vec3::Create(0) );
			TMEMSET( &m_tmpRender2[0], m_tmpRender2.size(), Vec3::Create(0) );
			
			// first do excessive rasterization
			RasterizeInstance( mi, 1.0f );
			
			// then do proper rasterization
			RasterizeInstance( mi, 0.0f );
		}
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
	case LTR_WT_PREXFORM:
		info->stage = "transforming spatial data";
		info->item_count = scene->m_meshInstances.size();
		break;
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

LTRBOOL ltr_MeshAddPart( ltr_Mesh* mesh, ltr_MeshPartInfo* mpinfo )
{
	ltr_MeshPart mp =
	{
		mpinfo->vertex_count,
		(u32) mesh->m_vpos.size(),
		mpinfo->index_count,
		(u32) mesh->m_indices.size(),
		mpinfo->tristrip
	};
	
	if( mpinfo->index_count < 3 && ( mpinfo->tristrip || mpinfo->index_count % 3 != 0 ) )
		return 0;
	
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
	return 1;
}

LTRBOOL ltr_MeshAddInstance( ltr_Mesh* mesh, ltr_MeshInstanceInfo* mii )
{
	ltr_MeshInstance* mi = new ltr_MeshInstance;
	mi->mesh = mesh;
	memcpy( mi->matrix.a, mii->matrix, sizeof(Mat4) );
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


