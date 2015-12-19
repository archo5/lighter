

#include "lighter_int.hpp"



LTRBOOL ltr_DefaultSizeFunc
(
	ltr_Config* config,
	const char* mesh_ident,
	size_t mesh_ident_size,
	const char* inst_ident,
	size_t inst_ident_size,
	float computed_surface_area,
	float inst_importance,
	u32 out_size[2]
)
{
	float size = config->global_size_factor * sqrtf( computed_surface_area ) * inst_importance;
	if( size < 1 )
		size = 1;
	u32 ui_size = ltr_NextPowerOfTwo( (u32) size );
	if( ui_size > config->max_lightmap_size )
		return 0;
	out_size[0] = ui_size;
	out_size[1] = ui_size;
	return 1;
}



void RasterizeInstance( ltr_MeshInstance* mi, float margin, Vec3* R1, Vec3* R2, Vec4* R3 )
{
	ltr_Mesh* mesh = mi->mesh;
	for( u32 part = 0; part < mesh->m_parts.size(); ++part )
	{
		const ltr_MeshPart& mp = mesh->m_parts[ part ];
		u32 O = mp.m_vertexOffset;
		const u32* indexBase = &mesh->m_indices[ mp.m_indexOffset ];
		
		for( u32 tri = 0; tri < mp.m_indexCount; tri += 3 )
		{
			u32 tridx1 = tri, tridx2 = tri + 1, tridx3 = tri + 2;
			
			tridx1 = O + indexBase[ tridx1 ];
			tridx2 = O + indexBase[ tridx2 ];
			tridx3 = O + indexBase[ tridx3 ];
			
			float samplearea = CalculateSampleArea(
				mi->m_ltex[ tridx1 ], mi->m_ltex[ tridx2 ], mi->m_ltex[ tridx3 ],
				mi->m_vpos[ tridx1 ], mi->m_vpos[ tridx2 ], mi->m_vpos[ tridx3 ]
			);
			
			RasterizeTriangle2D_x2_ex
			(
				R1, R2, R3, mi->lm_width, mi->lm_height, margin,
				mi->m_ltex[ tridx1 ], mi->m_ltex[ tridx2 ], mi->m_ltex[ tridx3 ],
				mi->m_vpos[ tridx1 ], mi->m_vpos[ tridx2 ], mi->m_vpos[ tridx3 ],
				mi->m_vnrm[ tridx1 ], mi->m_vnrm[ tridx2 ], mi->m_vnrm[ tridx3 ],
				V4( mi->m_vtex[ tridx1 ].x, mi->m_vtex[ tridx1 ].y, (float) part, samplearea ),
				V4( mi->m_vtex[ tridx2 ].x, mi->m_vtex[ tridx2 ].y, (float) part, samplearea ),
				V4( mi->m_vtex[ tridx3 ].x, mi->m_vtex[ tridx3 ].y, (float) part, samplearea )
			);
		}
	}
}


struct QueryIDAdder
{
	std::vector< int32_t >* dest;
	void operator () ( int32_t* ids, int32_t count )
	{
		dest->reserve( dest->size() + count );
		for( int32_t i = 0; i < count; ++i )
			dest->push_back( ids[ i ] );
	}
};

void ltr_Light::QueryMeshInsts( AABBTree& tree, std::vector< int32_t >& out )
{
	QueryIDAdder query = { &out };
	if( type == LTR_LT_POINT )
	{
		Vec3 range3 = V3( range );
		tree.Query( position - range3, position + range3, query );
	}
	else if( type == LTR_LT_SPOT )
	{
		// optimize?
		Vec3 range3 = V3( range );
		tree.Query( position - range3, position + range3, query );
	}
	else if( type == LTR_LT_DIRECT )
	{
		tree.GetAll( query );
	}
}


ltr_Scene* ltr_CreateScene()
{
	return new ltr_Scene;
}

void ltr_DestroyScene( ltr_Scene* scene )
{
	delete scene;
}


struct SceneAnyHitRayQuery : BaseRayQuery
{
	SceneAnyHitRayQuery( ltr_Scene* s, const Vec3& r0, const Vec3& r1 ) : hit(false), S(s), ray_end(r1)
	{
		SetRay( r0, r1 );
	}
	bool operator () ( int32_t* ids, int32_t count )
	{
		for( int32_t i = 0; i < count; ++i )
		{
			ltr_MeshInstance* mi = S->m_meshInstances[ ids[ i ] ];
			if( mi->m_shadow )
			{
				hit = mi->m_triTree.IntersectRay( ray_origin, ray_end );
				if( hit )
					return false;
			}
		}
		return true;
	}
	
	bool hit;
	ltr_Scene* S;
	Vec3 ray_end;
};

bool ltr_Scene::VisibilityTest( const Vec3& A, const Vec3& B )
{
	Vec3 diffnorm = ( B - A ).Normalized();
	Vec3 mA = A + diffnorm * SMALL_FLOAT;
	Vec3 mB = B - diffnorm * SMALL_FLOAT;
	
	SceneAnyHitRayQuery query( this, mA, mB );
	m_instTree.RayQuery( query );
	return query.hit;
}


struct SceneDistanceQuery
{
	SceneDistanceQuery( ltr_Scene* s, const Vec3& p ) : pos( p ), dist( MAX_PENUMBRA_SIZE ), S( s ){ RecalcBB(); }
	
	void operator () ( int32_t* ids, int32_t count )
	{
		for( int32_t i = 0; i < count; ++i )
		{
			ltr_MeshInstance* mi = S->m_meshInstances[ ids[ i ] ];
			if( mi->m_shadow )
			{
				float ndst = mi->m_triTree.GetDistance( pos, dist );
				if( ndst < dist )
				{
					dist = ndst;
					RecalcBB();
				}
			}
		}
	}
	
	FORCEINLINE void RecalcBB()
	{
		bbmin = pos - V3(dist);
		bbmax = pos + V3(dist);
	}
	
	Vec3 bbmin, bbmax;
	Vec3 pos;
	float dist;
	ltr_Scene* S;
};

float ltr_Scene::Distance( const Vec3& p )
{
	SceneDistanceQuery query( this, p );
	m_instTree.DynBBQuery( query );
	return query.dist;
}

float ltr_Scene::CalcInvShadowFactor( const Vec3& from, const Vec3& to, float k )
{
	Vec3 rd = ( to - from ).Normalized();
	float mint = 0.001f;
	float maxt = ( to - from ).Length();
	
	float res = 1.0f;
	for( float t = mint; t < maxt; )
	{
		float h = Distance( from + rd * t );
		if( h < 0.001f )
			return 0.0f;
		res = TMIN( res, h / TMIN( t * k, MAX_PENUMBRA_SIZE ) );
		h = TMIN( h, MAX_PENUMBRA_STEP );
		t += h;
	}
	return res;
}


struct SceneClosestHitRayQueryBSP : BaseRayQuery
{
	SceneClosestHitRayQueryBSP( ltr_Scene* s, const Vec3& r0, const Vec3& r1, Vec3* on ) : closest(2), outnormal(on), S(s), ray_end(r1)
	{
		SetRay( r0, r1 );
	}
	bool operator () ( int32_t* ids, int32_t count )
	{
		for( int32_t i = 0; i < count; ++i )
		{
			ltr_MeshInstance* mi = S->m_meshInstances[ ids[ i ] ];
			if( mi->m_shadow )
			{
				Vec3 nrm;
				float dist = mi->m_bspTree.IntersectRay( ray_origin, ray_end, &nrm );
				if( dist < closest )
				{
					closest = dist;
					if( outnormal )
						*outnormal = nrm;
				}
			}
		}
		return true;
	}
	
	float closest;
	Vec3* outnormal;
	ltr_Scene* S;
	Vec3 ray_end;
};

float ltr_Scene::DistanceTest( const Vec3& A, const Vec3& B, Vec3* outnormal )
{
	Vec3 diffnorm = ( B - A ).Normalized();
	Vec3 mA = A + diffnorm * SMALL_FLOAT;
	Vec3 mB = B - diffnorm * SMALL_FLOAT;
	
	SceneClosestHitRayQueryBSP query( this, mA, mB, outnormal );
	m_instTree.RayQuery( query );
	return query.closest;
}


struct SceneClosestHitRayQueryBBT : BaseRayQuery
{
	SceneClosestHitRayQueryBBT( ltr_Scene* s, const Vec3& r0, const Vec3& r1 ) : closest(2), S(s), ray_end(r1)
	{
		SetRay( r0, r1 );
	}
	bool operator () ( int32_t* ids, int32_t count )
	{
		for( int32_t i = 0; i < count; ++i )
		{
			ltr_MeshInstance* mi = S->m_meshInstances[ ids[ i ] ];
			if( mi->m_shadow )
			{
				float dist = mi->m_triTree.IntersectRayDist( ray_origin, ray_end, NULL );
				if( dist < closest )
					closest = dist;
			}
		}
		return true;
	}
	
	float closest;
	ltr_Scene* S;
	Vec3 ray_end;
};

float ltr_Scene::DistanceTestBBT( const Vec3& A, const Vec3& B )
{
	Vec3 diffnorm = ( B - A ).Normalized();
	Vec3 mA = A + diffnorm * SMALL_FLOAT;
	Vec3 mB = B - diffnorm * SMALL_FLOAT;
	
	SceneClosestHitRayQueryBBT query( this, mA, mB );
	m_instTree.RayQuery( query );
	return query.closest;
}


void ltr_Scene::Job_PreXForm_Inner( ltr_MeshInstance* mi )
{
	if( mi->m_samplecont )
		return;
	ltr_Mesh* mesh = mi->mesh;
	
	mi->m_vpos.resize( mesh->m_vpos.size() );
	mi->m_vnrm.resize( mesh->m_vnrm.size() );
	mi->m_vtex.resize( mesh->m_vtex1.size() );
	mi->m_ltex.resize( mesh->m_vtex2.size() );
	
	TransformPositions( VDATA( mi->m_vpos ), VDATA( mesh->m_vpos ), mesh->m_vpos.size(), mi->matrix );
	TransformNormals( VDATA( mi->m_vnrm ), VDATA( mesh->m_vnrm ), mesh->m_vnrm.size(), mi->matrix );
	
	// compute total area
	float total_area = 0.0f;
	for( u32 part = 0; part < mesh->m_parts.size(); ++part )
	{
		const ltr_MeshPart& mp = mesh->m_parts[ part ];
		const Vec3* vertexBase = &mi->m_vpos[ mp.m_vertexOffset ];
		const u32* indexBase = &mesh->m_indices[ mp.m_indexOffset ];
		
		for( u32 tri = 0; tri < mp.m_indexCount; tri += 3 )
		{
			u32 tridx1 = tri, tridx2 = tri + 1, tridx3 = tri + 2;
			tridx1 = indexBase[ tridx1 ];
			tridx2 = indexBase[ tridx2 ];
			tridx3 = indexBase[ tridx3 ];
			total_area += TriangleArea( vertexBase[ tridx1 ], vertexBase[ tridx2 ], vertexBase[ tridx3 ] );
		}
	}
	
	// call lightmap resizer
	u32 out_size[2] = { config.default_width, config.default_height };
	if( !config.size_fn( &config, mesh->m_ident.c_str(), mesh->m_ident.size(), mi->m_ident.c_str(), mi->m_ident.size(), total_area, mi->m_importance, out_size ) )
	{
		out_size[0] = config.default_width;
		out_size[1] = config.default_height;
	}
	mi->lm_width = out_size[0];
	mi->lm_height = out_size[1];
	
	// transform texcoords
	Vec2 lsize = Vec2::Create( (float) mi->lm_width, (float) mi->lm_height );
	for( size_t i = 0; i < mi->m_ltex.size(); ++i )
	{
		mi->m_vtex[ i ] = mesh->m_vtex1[ i ];
		mi->m_ltex[ i ] = mesh->m_vtex2[ i ] * lsize - Vec2::Create(0.5f);
	}
}

void ltr_Scene::Job_PreXForm( LTRWorker::IO* io )
{
	ltr_Scene* S = (ltr_Scene*) io->shared;
	S->Job_PreXForm_Inner( S->m_meshInstances[ io->i ] );
}

void ltr_Scene::Job_ColInfo_Inner( ltr_MeshInstance* mi )
{
	if( mi->m_samplecont )
	{
		mi->m_triTree.SetTris( NULL, 0 );
		return;
	}

	std::vector< Triangle > tris;
	
	ltr_Mesh* mesh = mi->mesh;
	
	for( u32 part = 0; part < mesh->m_parts.size(); ++part )
	{
		const ltr_MeshPart& mp = mesh->m_parts[ part ];
		if( !mp.m_shadow )
			continue;
		const Vec3* vertexBase = &mi->m_vpos[ mp.m_vertexOffset ];
		const u32* indexBase = &mesh->m_indices[ mp.m_indexOffset ];
		
		for( u32 tri = 0; tri < mp.m_indexCount; tri += 3 )
		{
			u32 tridx1 = tri, tridx2 = tri + 1, tridx3 = tri + 2;
			tridx1 = indexBase[ tridx1 ];
			tridx2 = indexBase[ tridx2 ];
			tridx3 = indexBase[ tridx3 ];
			Triangle T = { vertexBase[ tridx1 ], vertexBase[ tridx2 ], vertexBase[ tridx3 ] };
			if( T.CheckIsUseful() )
				tris.push_back( T );
		}
	}
	
	if( tris.size() )
		mi->m_bspTree.SetTriangles( VDATA( tris ), tris.size() );
	mi->m_triTree.SetTris( VDATA( tris ), tris.size() );
}

void ltr_Scene::Job_ColInfo( LTRWorker::IO* io )
{
	ltr_Scene* S = (ltr_Scene*) io->shared;
	S->Job_ColInfo_Inner( S->m_meshInstances[ io->i ] );
}

void ltr_Scene::Job_Samples_Inner( ltr_Scene* S, ltr_MeshInstance* mi )
{
	float corr_min_dot = cosf( config.max_correct_angle / 180.0f * (float) M_PI );
	if( mi->m_samplecont )
	{
		for( size_t i = 0; i < mi->m_lightmap.size(); ++i )
			mi->m_lightmap[ i ] = Vec3::CreateFromPtr( config.ambient_color );
		return;
	}
	
	Vec3Vector R1;
	Vec3Vector R2;
	Vec4Vector R3;
	
	R1.resize( mi->lm_width * mi->lm_height );
	R2.resize( mi->lm_width * mi->lm_height );
	R3.resize( mi->lm_width * mi->lm_height );
	
	TMEMSET( VDATA( R1 ), R1.size(), V3(0) );
	TMEMSET( VDATA( R2 ), R2.size(), V3(0) );
	TMEMSET( VDATA( R3 ), R3.size(), V4(0) );
	
	// first do excessive rasterization
//	RasterizeInstance( mi, 2.0f + SMALL_FLOAT, VDATA( R1 ), VDATA( R2 ), VDATA( R3 ) );
//	RasterizeInstance( mi, 1.5f + SMALL_FLOAT, VDATA( R1 ), VDATA( R2 ), VDATA( R3 ) );
//	RasterizeInstance( mi, 1.0f + SMALL_FLOAT, VDATA( R1 ), VDATA( R2 ), VDATA( R3 ) );
	RasterizeInstance( mi, 0.5f + SMALL_FLOAT, VDATA( R1 ), VDATA( R2 ), VDATA( R3 ) );
	
	// then do proper rasterization
	RasterizeInstance( mi, 0.0f, VDATA( R1 ), VDATA( R2 ), VDATA( R3 ) );
	
	// compress samples to packed array
	u32 sample_count = 0;
	for( size_t i = 0; i < R2.size(); ++i )
	{
		Vec3& N = R2[ i ];
		if( !N.IsZero() )
			sample_count++;
	}
	mi->m_samples_pos.reserve( sample_count );
	mi->m_samples_nrm.reserve( sample_count );
	mi->m_samples_loc.reserve( sample_count );
	mi->m_samples_radinfo.reserve( sample_count );
	for( size_t i = 0; i < R2.size(); ++i )
	{
		Vec3 N = R2[ i ];
		if( !N.IsZero() )
		{
			N = N.Normalized();
			Vec3 P = R1[ i ];
			Vec4 XD = R3[ i ];
			
			// move sample out of concave edges
			for( size_t mid = 0; mid < S->m_meshInstances.size(); ++mid )
				S->m_meshInstances[ mid ]->m_triTree.OffsetSample( P, N, sqrtf( XD.w ) );
			
			// move sample in front of overlapping faces
			if( config.max_correct_dist )
			{
				int itsleft = 100;
				float md = config.max_correct_dist;
				Vec3 PEnd = P + N * md;
				while( md > SMALL_FLOAT && itsleft --> 0 )
				{
					Vec3 hitnrm = {0,0,0};
					float q = DistanceTest( P, PEnd, &hitnrm );
					if( -Vec3Dot( hitnrm, N ) < corr_min_dot )
						break;
					if( q < SMALL_FLOAT )
						q = SMALL_FLOAT;
					P = P * ( 1.0f - q ) + PEnd * q;
					md *= ( 1.0f - q );
				}
			}
			
			mi->m_samples_pos.push_back( P );
			mi->m_samples_nrm.push_back( N );
			mi->m_samples_loc.push_back( i );
			mi->m_samples_radinfo.push_back( XD );
		}
	}
	mi->m_lightmap.resize( sample_count );
	for( size_t i = 0; i < mi->m_lightmap.size(); ++i )
		mi->m_lightmap[ i ] = Vec3::CreateFromPtr( config.ambient_color );
}

void ltr_Scene::Job_Samples( LTRWorker::IO* io )
{
	ltr_Scene* S = (ltr_Scene*) io->shared;
	S->Job_Samples_Inner( S, S->m_meshInstances[ io->i ] );
}


void ltr_Scene::Job_LMRender_Point_Inner( size_t i, dw_lmrender_data* data )
{
	ltr_MeshInstance* mi = data->mi;
	ltr_Light& light = *data->light;
	
	// MAY BE THREADED
	// for( size_t i = 0; i < mi->m_lightmap.size(); ++i )
	{
		Vec3& SP = mi->m_samples_pos[ i ];
		Vec3& SN = mi->m_samples_nrm[ i ];
		Vec3 sample2light = light.position - SP;
		float dist = sample2light.Length();
		if( dist )
			sample2light /= dist;
		float f_dist = pow( 1 - TMIN( 1.0f, dist / light.range ), light.power );
		float f_ndotl = TMAX( 0.0f, Vec3Dot( sample2light, SN ) );
		if( f_dist * f_ndotl <= 0 )
			return; // continue;
		float f_vistest = CalcInvShadowFactor( SP + SN * SAMPLE_SHADOW_OFFSET, light.position, light.light_radius );
		mi->m_lightmap[ i ] += light.color_rgb * ( f_dist * f_ndotl * f_vistest );
		if( config.generate_normalmap_data )
		{
			float factor = f_dist * f_vistest * CalcBrightness( light.color_rgb );
			if( factor > 0 )
			{
				ltr_LightContribSample contrib = { sample2light * factor, (int32_t) i };
				(*data->contribs)[ i ] = contrib;
			}
		}
	}
}

void ltr_Scene::Job_LMRender_Spot_Inner( size_t i, dw_lmrender_data* data )
{
	ltr_MeshInstance* mi = data->mi;
	ltr_Light& light = *data->light;
	float angle_out_rad = data->angle_out_rad;
//	float angle_in_rad = data->angle_in_rad;
	float angle_diff = data->angle_diff;
	
	// MAY BE THREADED
	// for( size_t i = 0; i < mi->m_lightmap.size(); ++i )
	{
		Vec3& SP = mi->m_samples_pos[ i ];
		Vec3& SN = mi->m_samples_nrm[ i ];
		Vec3 sample2light = light.position - SP;
		float dist = sample2light.Length();
		if( dist )
			sample2light /= dist;
		float f_dist = pow( 1 - TMIN( 1.0f, dist / light.range ), light.power );
		float f_ndotl = TMAX( 0.0f, Vec3Dot( sample2light, SN ) );
		float angle = acosf( TMIN( 1.0f, Vec3Dot( sample2light, -light.direction ) ) );
		float f_dir = TMAX( 0.0f, TMIN( 1.0f, ( angle - angle_out_rad ) / angle_diff ) );
		f_dir = pow( f_dir, light.spot_curve );
		if( f_dist * f_ndotl * f_dir <= 0 )
			return; // continue;
		float f_vistest = CalcInvShadowFactor( SP + SN * SAMPLE_SHADOW_OFFSET, light.position, light.light_radius );
		mi->m_lightmap[ i ] += light.color_rgb * ( f_dist * f_ndotl * f_dir * f_vistest );
		if( config.generate_normalmap_data )
		{
			float factor = f_dist * f_dir * f_vistest * CalcBrightness( light.color_rgb );
			if( factor > 0 )
			{
				ltr_LightContribSample contrib = { sample2light * factor, (int32_t) i };
				(*data->contribs)[ i ] = contrib;
			}
		}
	}
}

void ltr_Scene::Job_LMRender_Direct_Inner( size_t i, dw_lmrender_data* data )
{
	ltr_MeshInstance* mi = data->mi;
	ltr_Light& light = *data->light;
	
	// MAY BE THREADED
	// for( size_t i = 0; i < mi->m_lightmap.size(); ++i )
	{
		Vec3& SP = mi->m_samples_pos[ i ];
		Vec3& SN = mi->m_samples_nrm[ i ];
		
#if 0
		float f_ndotl = 0.0f;
		float f_vistest = 0.0f;
		Vec3 ray_origin = SP;
		
		float randoff = randf();//( sin( SP.x ) + sin( SP.y ) + sin( SP.z ) ) / 6 + 0.5f;
		
		for( int s = 0; s < light.shadow_sample_count; ++s )
		{
			Vec3 adjdir = Vec3::CreateSpiralDirVector( -light.direction, randoff, s, light.shadow_sample_count );
			adjdir = ( adjdir + (-light.direction) * tan( ( light.light_radius - 0.5f ) * float(M_PI) * 0.999f ) ).Normalized();
			
			f_ndotl += TMAX( 0.0f, Vec3Dot( adjdir, SN ) );
			if( VisibilityTest( ray_origin, ray_origin + adjdir * light.range ) )
				f_vistest += 1.0f;
		}
		f_ndotl /= light.shadow_sample_count;
		f_vistest /= light.shadow_sample_count;
		f_vistest = 1.0f - f_vistest;
#else
		float f_ndotl = TMAX( 0.0f, Vec3Dot( light.direction, SN ) );
		float f_vistest = CalcInvShadowFactor( SP + SN * SAMPLE_SHADOW_OFFSET, SP + light.direction * light.range, light.light_radius );
#endif
		mi->m_lightmap[ i ] += light.color_rgb * ( f_ndotl * f_vistest );
		if( config.generate_normalmap_data )
		{
			float factor = f_vistest * CalcBrightness( light.color_rgb );
			if( factor > 0 )
			{
				ltr_LightContribSample contrib = { light.direction * factor, (int32_t) i };
				(*data->contribs)[ i ] = contrib;
			}
		}
	}
}

void ltr_Scene::Job_LMRender_Point( LTRWorker::IO* io )
{
	ltr_Scene* S = (ltr_Scene*) io->shared;
	dw_lmrender_data* data = (dw_lmrender_data*) io->item;
	S->Job_LMRender_Point_Inner( io->i, data );
}

void ltr_Scene::Job_LMRender_Spot( LTRWorker::IO* io )
{
	ltr_Scene* S = (ltr_Scene*) io->shared;
	dw_lmrender_data* data = (dw_lmrender_data*) io->item;
	S->Job_LMRender_Spot_Inner( io->i, data );
}

void ltr_Scene::Job_LMRender_Direct( LTRWorker::IO* io )
{
	ltr_Scene* S = (ltr_Scene*) io->shared;
	dw_lmrender_data* data = (dw_lmrender_data*) io->item;
	S->Job_LMRender_Direct_Inner( io->i, data );
}

void ltr_Scene::Int_LMRender( ltr_Light& light, ltr_MeshInstance* mi )
{
	std::vector<ltr_LightContribSample> contribs;
	if( config.generate_normalmap_data )
	{
		contribs.resize( mi->m_lightmap.size() );
		ltr_LightContribSample nullsample = { V3(0), -1 };
		TMEMSET( VDATA( contribs ), contribs.size(), nullsample );
	}
	
	if( light.type == LTR_LT_POINT )
	{
		dw_lmrender_data data = { &contribs, mi, &light };
		
		m_worker.DoWork( this, &data, 0, mi->m_lightmap.size(), ltr_Scene::Job_LMRender_Point );
	}
	if( light.type == LTR_LT_SPOT )
	{
		float angle_out_rad = light.spot_angle_out / 180.0f * (float) M_PI, angle_in_rad = light.spot_angle_in / 180.0f * (float) M_PI;
		if( angle_in_rad == angle_out_rad )
			angle_in_rad -= SMALL_FLOAT;
		float angle_diff = angle_in_rad - angle_out_rad;
		dw_lmrender_data data = { &contribs, mi, &light, angle_out_rad, angle_in_rad, angle_diff };
		
		m_worker.DoWork( this, &data, 0, mi->m_lightmap.size(), ltr_Scene::Job_LMRender_Spot );
	}
	else if( light.type == LTR_LT_DIRECT )
	{
		dw_lmrender_data data = { &contribs, mi, &light };
		
		m_worker.DoWork( this, &data, 0, mi->m_lightmap.size(), ltr_Scene::Job_LMRender_Direct );
	}
	
	// commit contributions
	if( config.generate_normalmap_data )
	{
		for( size_t i = 0; i < contribs.size(); ++i )
		{
			if( contribs[ i ].sid >= 0 )
				mi->m_contrib.push_back( contribs[ i ] );
		}
	}
}

void ltr_Scene::Int_RDGenLinks()
{
#ifdef LTRDEBUG
	double t0 = ltr_gettime();
#endif
	
	for( size_t mi_id = 0; mi_id < m_meshInstances.size(); ++mi_id )
	{
		ltr_MeshInstance* mi = m_meshInstances[ mi_id ];
		for( size_t i = 0; i < mi->m_lightmap.size(); ++i )
		{
			ltr_Mesh* mesh = mi->mesh;
			const Vec3& P = mi->m_samples_pos[ i ];
			const Vec3& N = mi->m_samples_nrm[ i ].Normalized();
			
			ltr_RadSampleGeom RSG = { P, N };
			m_radSampleGeoms.push_back( RSG );
			
			Vec3 s_diff = {1,1,1};
			Vec3 s_emit = mi->m_lightmap[ i ];
			float s_area = mi->m_samplecont ? 0 : mi->m_samples_radinfo[ i ].w;
			if( mi->m_samplecont )
			{
				s_diff = Vec3::Create(0);
			}
			else if( config.sample_fn )
			{
				const size_t L = mi->m_samples_loc[ i ];
				int LMX = L % mi->lm_width;
				int LMY = L / mi->lm_width;
				ltr_SampleRequest REQ =
				{
					{ P.x, P.y, P.z },
					{ N.x, N.y, N.z },
					mi->m_samples_radinfo[ i ].x, mi->m_samples_radinfo[ i ].y,
					( LMX + 0.5f ) / mi->lm_width, ( LMY + 0.5f ) / mi->lm_height,
					(uint32_t) mi->m_samples_radinfo[ i ].z, // should not be interpolated between triangles
					mesh->m_ident.c_str(), mesh->m_ident.size(),
					mi->m_ident.c_str(), mi->m_ident.size(),
					{ 1, 1, 1 },
					{ 0, 0, 0 },
				};
				if( config.sample_fn( &config, &REQ ) )
				{
					s_diff = Vec3::CreateFromPtr( REQ.out_diffuse_color );
					s_emit += Vec3::CreateFromPtr( REQ.out_emissive_color );
				}
			}
			
			ltr_RadSampleColors RSC = { s_diff, s_emit, s_emit, Vec3::Create(0), s_area };
			m_radSampleColors.push_back( RSC );
		}
	}
	
#ifdef LTRDEBUG
	double t1 = ltr_gettime();
	printf( "RAD sample concat: %g seconds\n", t1 - t0 );
	
	double t2 = ltr_gettime();
#endif
	
	for( size_t i = 0; i < m_radSampleGeoms.size(); ++i )
	{
		m_radLinkMap.push_back( m_radLinks.size() ); // offset
		m_radLinkMap.push_back( 0 ); // count
		
		const Vec3& A_SP = m_radSampleGeoms[ i ].pos;
		const Vec3& A_SN = m_radSampleGeoms[ i ].normal;
		for( size_t j = i + 1; j < m_radSampleGeoms.size(); ++j )
		{
			const Vec3& B_SP = m_radSampleGeoms[ j ].pos;
			const Vec3& B_SN = m_radSampleGeoms[ j ].normal;
			
			Vec3 posdiff = B_SP - A_SP;
			float dotA = Vec3Dot( A_SN, posdiff );
			float dotB = Vec3Dot( B_SN, -posdiff );
			if( dotA <= SMALL_FLOAT || dotB <= SMALL_FLOAT )
				continue;
			
			float lensq = posdiff.LengthSq();
			float factor = dotA * dotB / ( lensq * lensq * (float) M_PI );
			if( factor < SMALL_FLOAT )
				continue;
			
			if( VisibilityTest( A_SP, B_SP ) == false )
				continue;
			
			// add sample
			ltr_RadLink RL = { (uint32_t) j, factor };
			m_radLinks.push_back( RL );
			m_radLinkMap.back()++;
		}
	}
	
#ifdef LTRDEBUG
	double t3 = ltr_gettime();
	printf( "RAD sample processing: %g seconds\n", t3 - t2 );
#endif
}

void ltr_Scene::Int_RDBounce()
{
	for( size_t i = 0; i < m_radLinkMap.size() / 2; ++i )
	{
		size_t sidbegin = m_radLinkMap[ i*2+0 ];
		size_t sidend = sidbegin + m_radLinkMap[ i*2+1 ];
		for( size_t sid = sidbegin; sid < sidend; ++sid )
		{
			size_t j = m_radLinks[ sid ].other;
			float f = m_radLinks[ sid ].factor;
			
			ltr_RadSampleColors& A_RSC = m_radSampleColors[ i ];
			ltr_RadSampleColors& B_RSC = m_radSampleColors[ j ];
			
			Vec3 A_out = ( A_RSC.outputEnergy * A_RSC.diffuseColor ) * A_RSC.area * f;
			Vec3 B_out = ( B_RSC.outputEnergy * B_RSC.diffuseColor ) * B_RSC.area * f;
			
			A_RSC.totalLight += B_out;
			B_RSC.totalLight += A_out;
			
			A_RSC.inputEnergy += B_out;
			B_RSC.inputEnergy += A_out;
		}
	}
	for( size_t i = 0; i < m_radSampleColors.size(); ++i )
	{
		ltr_RadSampleColors& RSC = m_radSampleColors[ i ];
		RSC.outputEnergy = RSC.inputEnergy;
		RSC.inputEnergy = Vec3::Create(0);
	}
}

void ltr_Scene::Job_AORender_Inner( ltr_MeshInstance* mi, size_t i )
{
//	float ao_divergence = config.ao_divergence * 0.5f + 0.5f;
	float ao_distance = config.ao_distance,
		ao_falloff = config.ao_falloff,
		ao_multiplier = config.ao_multiplier,
		ao_effect = config.ao_effect;
	int num_samples = config.ao_num_samples;
	Vec3 ao_color = Vec3::CreateFromPtr( config.ao_color_rgb );
	
	// MAY BE THREADED
	// for( size_t i = 0; i < mi->m_lightmap.size(); ++i )
	{
		Vec3& SP = mi->m_samples_pos[ i ];
		Vec3& SN = mi->m_samples_nrm[ i ];
		Vec3& OutColor = mi->m_lightmap[ i ];
		
		Vec3 ray_origin = SP + SN * ( SMALL_FLOAT * 2 );
		
		float ao_factor = 0;
		float randoff = randf();
		for( int s = 0; s < num_samples; ++s )
		{
			Vec3 ray_dir = Vec3::CreateSpiralDirVector( SN, randoff, s, num_samples ) * ao_distance;
			float hit = DistanceTestBBT( ray_origin, ray_origin + ray_dir );
			if( hit < 1.0f )
				ao_factor += 1.0f - hit;
		}
		ao_factor /= num_samples;
		ao_factor = TMIN( ao_factor * ao_multiplier, 1.0f );
		if( ao_falloff )
			ao_factor = pow( ao_factor, ao_falloff );
		
		if( ao_effect >= 0 )
		{
			OutColor = OutColor * ( 1 - ao_factor * ( 1 - ao_effect ) ) + ao_color * ao_factor;
		}
		else
		{
			OutColor = TLERP( OutColor, TLERP( ao_color, ao_color * OutColor, -ao_effect ), ao_factor );
		}
	}
}

void ltr_Scene::Job_AORender( LTRWorker::IO* io )
{
	ltr_Scene* S = (ltr_Scene*) io->shared;
	S->Job_AORender_Inner( (ltr_MeshInstance*) io->item, io->i );
}

void ltr_Scene::Int_Finalize()
{
	for( size_t mid = 0; mid < m_meshInstances.size(); ++mid )
	{
		ltr_MeshInstance* mi = m_meshInstances[ mid ];
		if( mi->m_samplecont )
		{
			for( size_t i = 0; i < mi->m_lightmap.size(); ++i )
			{
				m_samples[ i ].out_color[0] = mi->m_lightmap[ i ].x;
				m_samples[ i ].out_color[1] = mi->m_lightmap[ i ].y;
				m_samples[ i ].out_color[2] = mi->m_lightmap[ i ].z;
			}
			continue;
		}
		ltr_Mesh* mesh = mi->mesh;
		
		float* image_rgb = new float[ mi->lm_width * mi->lm_height * 3 ];
		char* image_set = new char[ mi->lm_width * mi->lm_height * 2 ];
		TMEMSET( image_rgb, mi->lm_width * mi->lm_height * 3, 0.0f );
		TMEMSET( image_set, mi->lm_width * mi->lm_height * 2, (char) 0 );
		for( size_t i = 0; i < mi->m_samples_loc.size(); ++i )
		{
			size_t loc = mi->m_samples_loc[ i ] * 3;
			Vec3& LMC = mi->m_lightmap[ i ];
			image_rgb[ loc+0 ] = LMC.x;
			image_rgb[ loc+1 ] = LMC.y;
			image_rgb[ loc+2 ] = LMC.z;
			image_set[ mi->m_samples_loc[ i ] ] = 1;
		}
		
		// extrapolate for interpolation
		u32 w = mi->lm_width, h = mi->lm_height;
#define PX_ISSET( x, y ) ( (x) < w && (y) < h && image_set[ (x) + (y) * w ] )
#define PX_GETCOL( x, y ) Vec3::Create(image_rgb[((x)+(y)*w)*3+0],image_rgb[((x)+(y)*w)*3+1],image_rgb[((x)+(y)*w)*3+2])
		for( int its = 1; its <= 3; ++its )
		{
			memcpy( image_set + w*h, image_set, w*h );
			for( u32 y = 0; y < h; ++y )
			{
				for( u32 x = 0; x < w; ++x )
				{
					if( image_set[ x + y * w ] )
						continue;
					
					Vec3 col = {0,0,0};
					int count = 0;
					// H
					if( PX_ISSET( x-1, y ) )
					{
						if( PX_ISSET( x-2, y ) )
							col += TLERP( PX_GETCOL( x-2, y ), PX_GETCOL( x-1, y ), 2.0f );
						else
							col += PX_GETCOL( x-1, y );
						count++;
					}
					if( PX_ISSET( x+1, y ) )
					{
						if( PX_ISSET( x+2, y ) )
							col += TLERP( PX_GETCOL( x+2, y ), PX_GETCOL( x+1, y ), 2.0f );
						else
							col += PX_GETCOL( x+1, y );
						count++;
					}
					// V
					if( PX_ISSET( x, y-1 ) )
					{
						if( PX_ISSET( x, y-2 ) )
							col += TLERP( PX_GETCOL( x, y-2 ), PX_GETCOL( x, y-1 ), 2.0f );
						else
							col += PX_GETCOL( x, y-1 );
						count++;
					}
					if( PX_ISSET( x, y+1 ) )
					{
						if( PX_ISSET( x, y+2 ) )
							col += TLERP( PX_GETCOL( x, y+2 ), PX_GETCOL( x, y+1 ), 2.0f );
						else
							col += PX_GETCOL( x, y+1 );
						count++;
					}
					
					if( count )
					{
						col /= (float) count;
						u32 off = ( x + y * w ) * 3;
						image_rgb[ off+0 ] = TMAX( col.x, 0.0f );
						image_rgb[ off+1 ] = TMAX( col.y, 0.0f );
						image_rgb[ off+2 ] = TMAX( col.z, 0.0f );
						image_set[ x + y * w + w*h ] = 1;
					}
				}
			}
			memcpy( image_set, image_set + w*h, w*h );
		}
#undef PX_ISSET
#undef PX_GETCOL
		
		if( config.blur_size )
		{
			int blur_ext = (int) ceil( config.blur_size );
			int blurbuf_size = ( TMAX( mi->lm_width, mi->lm_height ) + blur_ext * 2 ) * 3;
			int kernel_size = blur_ext * 2 + 1;
			float* image_tmp_rgb = new float[ mi->lm_width * mi->lm_height * 3 + blurbuf_size + kernel_size ];
			float* blur_tmp_rgb = image_tmp_rgb + mi->lm_width * mi->lm_height * 3;
			float* kernel_rgb = blur_tmp_rgb + blurbuf_size;
			
			Generate_Gaussian_Kernel( kernel_rgb, blur_ext, config.blur_size );
			Convolve_Transpose( image_rgb, image_tmp_rgb, mi->lm_width, mi->lm_height, blur_ext, kernel_rgb, blur_tmp_rgb );
			Convolve_Transpose( image_tmp_rgb, image_rgb, mi->lm_height, mi->lm_width, blur_ext, kernel_rgb, blur_tmp_rgb );
			
			delete [] image_tmp_rgb;
		}
		
		if( config.ds2x )
		{
			int hw = TMAX( mi->lm_width / 2, u32(1) );
			int hh = TMAX( mi->lm_height / 2, u32(1) );
			float* image_ds2x = new float[ hw * hh * 3 ];
			
			Downsample2X( image_ds2x, hw, hh, image_rgb, mi->lm_width, mi->lm_height );
			
			delete [] image_rgb;
			image_rgb = image_ds2x;
			mi->lm_width = hw;
			mi->lm_height = hh;
		}
		
		float* normals_xyzf = NULL;
		if( config.generate_normalmap_data )
		{
			std::vector<ltr_TmpContribSum> contribs;
			
			// - prepare temp. data store
			contribs.resize( mi->m_samples_loc.size() );
			ltr_TmpContribSum tmpl = { V3(0), 1, 0 };
			TMEMSET( VDATA( contribs ), contribs.size(), tmpl );
			
			// - gather contributions and average direction
			for( size_t i = 0; i < mi->m_contrib.size(); ++i )
			{
				contribs[ mi->m_contrib[ i ].sid ].normal += mi->m_contrib[ i ].normal;
				contribs[ mi->m_contrib[ i ].sid ].count++;
			}
			for( size_t i = 0; i < contribs.size(); ++i )
			{
				if( contribs[ i ].count > 0 )
					contribs[ i ].normal /= contribs[ i ].count;
			}
			
			// - calculate focus
			for( size_t i = 0; i < mi->m_contrib.size(); ++i )
			{
				ltr_TmpContribSum& TCS = contribs[ mi->m_contrib[ i ].sid ];
				Vec3 cn = mi->m_contrib[ i ].normal;
				float dot = 1 - ( 1 - TMAX( Vec3Dot( TCS.normal.Normalized(), cn.Normalized() ), 0.0f ) ) / TMAX( TCS.normal.Length() / cn.Length(), 1.0f );
				TCS.mindot = TMIN( TCS.mindot, dot );
			}
			
			// - reorder data into a texture
			normals_xyzf = new float[ mi->lm_width * mi->lm_height * 4 ];
			TMEMSET( normals_xyzf, mi->lm_width * mi->lm_height * 4, 0.0f );
			for( size_t i = 0; i < mi->m_samples_loc.size(); ++i )
			{
				ltr_TmpContribSum& TCS = contribs[ i ];
				size_t loc = mi->m_samples_loc[ i ] * 4;
				Vec3 N = TCS.normal.Normalized();
				normals_xyzf[ loc + 0 ] = N.x;
				normals_xyzf[ loc + 1 ] = N.y;
				normals_xyzf[ loc + 2 ] = N.z;
				normals_xyzf[ loc + 3 ] = TCS.mindot;
			}
		}
		
		ltr_WorkOutput wo =
		{
			mid,
			mesh->m_ident.c_str(),
			mesh->m_ident.size(),
			mi->m_ident.c_str(),
			mi->m_ident.size(),
			image_rgb,
			normals_xyzf,
			mi->lm_width,
			mi->lm_height
		};
		
		m_workOutput.push_back( wo );
	}
}

void ltr_Scene::Job_MainProc( LTRWorker::IO* io )
{
	ltr_Scene* S = (ltr_Scene*) io->shared;
	
	if( S->m_meshInstances.size() )
	{
		S->m_workStage = "transforming spatial data";
		S->m_workCompletion = 0;
		S->m_worker.DoWork( S, NULL, 0, S->m_meshInstances.size(), ltr_Scene::Job_PreXForm );
		
		S->m_workStage = "generating data structures";
		S->m_workCompletion = 0;
		S->m_worker.DoWork( S, NULL, 0, S->m_meshInstances.size(), ltr_Scene::Job_ColInfo );
		
		S->m_workCompletion = 0.99f;
		std::vector< AABB3 > bbs;
		for( size_t i = 0; i < S->m_meshInstances.size(); ++i )
		{
	//		S->m_meshInstances[ i ]->m_triTree.m_bbTree.Dump();
			AABBTree::Node& N = S->m_meshInstances[ i ]->m_triTree.m_bbTree.m_nodes[ 0 ];
			AABB3 bb = { N.bbmin, N.bbmax };
			bbs.push_back( bb );
		}
		S->m_instTree.SetAABBs( &bbs[0], bbs.size() );
		
		S->m_workStage = "generating samples";
		S->m_workCompletion = 0;
		S->m_worker.DoWork( S, NULL, 0, S->m_meshInstances.size(), ltr_Scene::Job_Samples );
		
		if( S->m_lights.size() )
		{
			S->m_workStage = "rendering lightmaps";
			S->m_workCompletion = 0;
			std::vector< int32_t > mesh_inst_ids;
			for( size_t lt_id = 0; lt_id < S->m_lights.size(); ++lt_id )
			{
				mesh_inst_ids.clear();
				mesh_inst_ids.push_back( 0 ); // sample container
				S->m_lights[ lt_id ].QueryMeshInsts( S->m_instTree, mesh_inst_ids );
				std::sort( mesh_inst_ids.begin(), mesh_inst_ids.end() );
				int32_t previd = -1;
				for( size_t i = 0; i < mesh_inst_ids.size(); ++i )
				{
					if( i > 0 && mesh_inst_ids[ i ] == previd )
						continue;
					previd = mesh_inst_ids[ i ];
					S->Int_LMRender( S->m_lights[ lt_id ], S->m_meshInstances[ previd ] );
				}
				S->m_workCompletion = ( lt_id + 1 ) / float(S->m_lights.size());
			}
		}
		
		if( S->config.bounce_count )
		{
			S->m_workStage = "calculating radiosity";
			S->m_workCompletion = 0;
			S->Int_RDGenLinks();
			
			S->m_workStage = "bouncing light";
			S->m_workCompletion = 0;
			for( int b = 0; b < S->config.bounce_count; ++b )
			{
				S->Int_RDBounce();
				S->m_workCompletion = ( b + 1 ) / float(S->config.bounce_count);
			}
			
			S->m_workStage = "committing radiosity";
			S->m_workCompletion = 0;
			size_t sp = 0;
			for( size_t mi_id = 0; mi_id < S->m_meshInstances.size(); ++mi_id )
			{
				ltr_MeshInstance* mi = S->m_meshInstances[ mi_id ];
				for( size_t i = 0; i < mi->m_lightmap.size(); ++i )
				{
					mi->m_lightmap[ i ] = S->m_radSampleColors[ sp++ ].totalLight;
				}
				S->m_workCompletion = ( mi_id + 1 ) / float(S->m_meshInstances.size());
			}
		}
		
		if( S->config.ao_distance )
		{
			S->m_workStage = "rendering ambient occlusion";
			S->m_workCompletion = 0;
			for( size_t mi_id = 0; mi_id < S->m_meshInstances.size(); ++mi_id )
			{
				ltr_MeshInstance* mi = S->m_meshInstances[ mi_id ];
				S->m_worker.DoWork( S, mi, 0, mi->m_lightmap.size(), ltr_Scene::Job_AORender );
				S->m_workCompletion = ( mi_id + 1 ) / float(S->m_meshInstances.size());
			}
		}
		
		S->m_workStage = "exporting lightmaps";
		S->m_workCompletion = 0;
		S->Int_Finalize();
	}
	
	S->m_workStage = NULL;
	S->m_workCompletion = 1;
}

void ltr_Start( ltr_Scene* scene )
{
	scene->m_worker.Init( TMAX( TMIN( scene->config.max_num_threads, scene->m_num_cpus ), 1 ) );
	scene->m_rootWorker.Init( 1 );
	scene->m_rootWorker.DoWork( scene, NULL, 0, 1, ltr_Scene::Job_MainProc, false );
}

void ltr_Abort( ltr_Scene* scene )
{
	// TODO
}

LTRBOOL ltr_GetStatus( ltr_Scene* scene, ltr_WorkStatus* wsout )
{
	wsout->completion = scene->m_workCompletion;
	wsout->stage = scene->m_workStage ? scene->m_workStage : "finished";
	return scene->m_workStage != NULL;
}

void ltr_Sleep( int ms )
{
	ltrthread_sleep( (uint32_t) ms );
}

void ltr_GetConfig( ltr_Config* cfg, ltr_Scene* opt_scene )
{
	if( opt_scene )
	{
		memcpy( cfg, &opt_scene->config, sizeof(*cfg) );
		return;
	}
	
	cfg->userdata = NULL;
	cfg->size_fn = ltr_DefaultSizeFunc;
	
	cfg->max_num_threads = 0x7fff;
	cfg->max_tree_memory = 128 * 1024 * 1024;
	cfg->max_lightmap_size = 1024;
	cfg->default_width = 64;
	cfg->default_height = 64;
	cfg->global_size_factor = 4;
	cfg->max_correct_dist = 0.1f;
	cfg->max_correct_angle = 60;
	
	LTR_VEC3_SET( cfg->clear_color, 0, 0, 0 );
	LTR_VEC3_SET( cfg->ambient_color, 0, 0, 0 );
	
	cfg->bounce_count = 0;
	cfg->sample_fn = NULL;
	
	cfg->ao_distance = 0;
	cfg->ao_multiplier = 1.2f;
	cfg->ao_falloff = 1;
	cfg->ao_effect = 0;
	cfg->ao_divergence = 0;
	LTR_VEC3_SET( cfg->ao_color_rgb, 0, 0, 0 );
	cfg->ao_num_samples = 17;
	
	cfg->blur_size = 0.5f;
	cfg->ds2x = 0;
	cfg->generate_normalmap_data = 0;
}

LTRCODE ltr_SetConfig( ltr_Scene* scene, ltr_Config* cfg )
{
	memcpy( &scene->config, cfg, sizeof(*cfg) );
	return 1;
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
		mpinfo->shadow
	};
	
	if( mpinfo->index_count < 3 && mpinfo->index_count % 3 != 0 )
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
	mi->m_importance = mii->importance;
	mi->m_shadow = !!mii->shadow;
	mi->m_samplecont = false;
	if( mii->ident )
		mi->m_ident.assign( mii->ident, mii->ident_size );
	memcpy( mi->matrix.a, mii->matrix, sizeof(Mat4) );
	mesh->m_scene->m_meshInstances.push_back( mi );
	mi->lm_width = 128;
	mi->lm_height = 128;
	return 1;
}

void ltr_LightAdd( ltr_Scene* scene, ltr_LightInfo* li )
{
	ltr_Light L =
	{
		li->type,
		Vec3::CreateFromPtr( li->position ),
		Vec3::CreateFromPtr( li->direction ).Normalized(),
		Vec3::CreateFromPtr( li->up_direction ).Normalized(),
		Vec3::CreateFromPtr( li->color_rgb ),
		li->range,
		li->power,
		li->light_radius,
		TMAX( 1, li->shadow_sample_count ),
		li->spot_angle_out,
		li->spot_angle_in,
		li->spot_curve,
	};
	L.samples.reserve( L.shadow_sample_count );
	float randoff = randf();
	for( int s = 0; s < L.shadow_sample_count; ++s )
	{
		if( L.type == LTR_LT_DIRECT )
		{
		//	Vec3 adjdir = Vec3::CreateRandomVectorDirDvg( -L.direction, L.light_radius );
			Vec3 adjdir = Vec3::CreateSpiralDirVector( -L.direction, randoff, s, L.shadow_sample_count );
			adjdir = ( adjdir + (-L.direction) * tan( ( L.light_radius - 0.5f ) * float(M_PI) * 0.999f ) ).Normalized();
			
			L.samples.push_back( adjdir );
		}
		else
		{
			// TODO
		}
	}
	scene->m_lights.push_back( L );
}

void ltr_SampleAdd( ltr_Scene* scene, ltr_SampleInfo* si )
{
	ltr_SampleInfo S = *si;
	S.out_color[0] = 0;
	S.out_color[1] = 0;
	S.out_color[2] = 0;
	scene->m_sampleMI.m_samples_pos.push_back( V3P( S.position ) );
	scene->m_sampleMI.m_samples_nrm.push_back( V3P( S.normal ) );
	scene->m_sampleMI.m_lightmap.push_back( V3(0) );
	scene->m_samples.push_back( S );
}


void ltr_GetWorkOutputInfo( ltr_Scene* scene, ltr_WorkOutputInfo* woutinfo )
{
	woutinfo->lightmap_count = scene->m_workOutput.size();
	woutinfo->sample_count = scene->m_samples.size();
	woutinfo->samples = VDATA( scene->m_samples );
}

LTRBOOL ltr_GetWorkOutput( ltr_Scene* scene, u32 which, ltr_WorkOutput* wout )
{
	if( which >= scene->m_workOutput.size() )
		return 0;
	*wout = scene->m_workOutput[ which ];
	return 1;
}


u32 ltr_NextPowerOfTwo( u32 x )
{
	x--;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return x + 1;
}


