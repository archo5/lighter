

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "lighter.h"
#include "lighter_int.hpp"


#define LTR_VERIFY( x ) printf( "VERIFY: %s : %d\n", #x, x )

template< typename T > T safe_div( T a, T b ){ return b == T(0) ? T(0) : a / b; }

static double ltr_curtime()
{
#ifdef __linux
	struct timespec ts;
	clock_gettime( CLOCK_MONOTONIC, &ts );
	return (double) ts.tv_sec + 0.000000001 * (double) ts.tv_nsec;
#else
	clock_t clk = clock();
	return (double)( clk ) / (double)( CLOCKS_PER_SEC );
#endif
}


struct testmesh
{
	Vec3Vector positions;
	Vec3Vector normals;
	Vec2Vector texcoords;
	U32Vector indices;
	
	u32 InsertVertex( const Vec3& P, const Vec3& N, const Vec2& T )
	{
		for( size_t i = 0; i < positions.size(); ++i )
		{
			if( positions[i] == P && normals[i] == N && texcoords[i] == T )
				return i;
		}
		positions.push_back( P );
		normals.push_back( N );
		texcoords.push_back( T );
		return positions.size() - 1;
	}
	
	bool Load( const char* file )
	{
		FILE* fh = fopen( file, "r" );
		if( !fh )
			return false;
		
		Vec3Vector plist;
		Vec3Vector nlist;
		Vec2Vector tlist;
		
		int i[9];
		Vec3 tmp3;
		Vec2 tmp2;
		char bfr[5] = {0};
		while( !feof( fh ) && fscanf( fh, "%4s", bfr ) )
		{
			if( !strcmp( bfr, "v" ) )
			{
				if( fscanf( fh, "%f %f %f\n", &tmp3.x, &tmp3.y, &tmp3.z ) )
					plist.push_back( tmp3 );
			}
			else if( !strcmp( bfr, "vt" ) )
			{
				if( fscanf( fh, "%f %f\n", &tmp2.x, &tmp2.y ) )
				{
					tmp2.y = 1 - tmp2.y;
					tlist.push_back( tmp2 );
				}
			}
			else if( !strcmp( bfr, "vn" ) )
			{
				if( fscanf( fh, "%f %f %f\n", &tmp3.x, &tmp3.y, &tmp3.z ) )
					nlist.push_back( tmp3 );
			}
			else if( !strcmp( bfr, "f" ) )
			{
				if( fscanf( fh, "%d/%d/%d %d/%d/%d %d/%d/%d\n", i+0, i+1, i+2, i+3, i+4, i+5, i+6, i+7, i+8 ) )
				{
					for( int j = 0; j < 9; ++j )
						i[j]--;
					Vec3 v1p = plist[i[0]], v2p = plist[i[3]], v3p = plist[i[6]];
					Vec2 v1t = tlist[i[1]], v2t = tlist[i[4]], v3t = tlist[i[7]];
					Vec3 v1n = nlist[i[2]], v2n = nlist[i[5]], v3n = nlist[i[8]];
					
					indices.push_back( InsertVertex( v1p, v1n, v1t ) );
					indices.push_back( InsertVertex( v2p, v2n, v2t ) );
					indices.push_back( InsertVertex( v3p, v3n, v3t ) );
				}
			}
			bfr[0] = 0;
		}
		
		printf( "loaded mesh '%s' - %d verts, %d tris (P %d N %d T %d)\n",
			file, (int) positions.size(), (int) indices.size() / 3,
			(int) plist.size(), (int) nlist.size(), (int) tlist.size() );
		
		fclose( fh );
		
		return true;
	}
};


void dumpimg( const char* file, float* img_rgb, u32 width, u32 height, float top = 1.0f )
{
	float mult = 255.0f / top;
	
	FILE* fh = fopen( file, "wb" );
	if( !fh )
	{
		printf( "FAILED TO WRITE TO FILE: %s\n", file );
		return;
	}
	unsigned char header[18]={0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	
    header[12] = width         & 0xFF;
    header[13] = ( width >> 8)  & 0xFF;
    header[14] = (height)       & 0xFF; 
    header[15] = (height >> 8)  & 0xFF;
    header[16] = 24;
    
    fwrite( header, 1, 18, fh );
    for( u32 y = height; y > 0; )
    {
    	--y;
    	for( u32 x = 0; x < width; ++x )
    	{
    		u32 off = ( x + width * y ) * 3;
    		unsigned char btwr[3] =
    		{
    			(unsigned char) TMIN( img_rgb[ off+2 ] * mult, 255.0f ),
    			(unsigned char) TMIN( img_rgb[ off+1 ] * mult, 255.0f ),
    			(unsigned char) TMIN( img_rgb[ off+0 ] * mult, 255.0f )
    		};
    		fwrite( btwr, 1, 3, fh );
    	}
    }
    fclose( fh );
}


void DOWORK( ltr_Scene* scene )
{
	ltr_WorkStatus wstatus;
	double ta = ltr_curtime();
	ltr_Start( scene );
	printf( "\n" );
	while( ltr_GetStatus( scene, &wstatus ) )
	{
		printf( "\r%s ... %d%%", wstatus.stage, int(wstatus.completion * 100) );
		ltr_Sleep( 50 );
	}
	double tb = ltr_curtime();
	printf( "\n\ntime taken: %f\n", (tb-ta) );
}


void testfunc_basic()
{
	puts( "LIGHTER test [basic]" );
	
	ltr_Scene* scene = ltr_CreateScene();
	
	// MESH 1
	ltr_Mesh* mesh1 = ltr_CreateMesh( scene, "mesh1", strlen("mesh1") );
	
	// MESH 1 PART 1
	float m1_p1_positions[] = { -1, -1, 0,  1, -1, 0,  1, 1, 0,  -1, 1, 0 };
	float m1_p1_normals[] = { 0, 0, 1,  0, 0, 1,  0, 0, 1,  0, 0, 1 };
	float m1_p1_texcoords[] = { 0.1f, 0.1f,  0.9f, 0.1f,  0.9f, 0.9f,  0.1f, 0.9f };
	u32 m1_p1_indices[] = { 0, 1, 2,  2, 3, 0 };
	ltr_MeshPartInfo m1part1 =
	{
		m1_p1_positions, m1_p1_normals, m1_p1_texcoords, m1_p1_texcoords,
		sizeof(float)*3, sizeof(float)*3, sizeof(float)*2, sizeof(float)*2,
		m1_p1_indices, 4, 6, 1
	};
	LTR_VERIFY( ltr_MeshAddPart( mesh1, &m1part1 ) );
	
	// MESH 1 INSTANCE 1
	ltr_MeshInstanceInfo m1inst1;
	memset( &m1inst1, 0, sizeof(m1inst1) );
	m1inst1.importance = 1;
	m1inst1.shadow = 1;
	float matrix1[16] = { 1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
	memcpy( m1inst1.matrix, matrix1, sizeof(matrix1) );
	LTR_VERIFY( ltr_MeshAddInstance( mesh1, &m1inst1 ) );
	
	// MESH 1 INSTANCE 2
	ltr_MeshInstanceInfo m1inst2;
	memset( &m1inst2, 0, sizeof(m1inst2) );
	m1inst2.importance = 1;
	m1inst2.shadow = 1;
	float c = cos( 0.5f ) * 0.2f, s = sin( 0.5f ) * 0.2f;
	float matrix2[16] = { c, s, 0, 0,  -s, c, 0, 0,  0, 0, 0.2f, 0,  0, 0, 0.5f, 1 };
	memcpy( m1inst2.matrix, matrix2, sizeof(matrix2) );
	LTR_VERIFY( ltr_MeshAddInstance( mesh1, &m1inst2 ) );
	
	// MESH 1 INSTANCE 3
	ltr_MeshInstanceInfo m1inst3;
	memset( &m1inst3, 0, sizeof(m1inst3) );
	m1inst3.importance = 1;
	m1inst3.shadow = 1;
	float matrix3[16] = { 1.9f, 0, 0, 0,  0, 1.9f, 0, 0,  0, 0, 1.9f, 0,  0, 0, -0.5f, 1 };
	memcpy( m1inst3.matrix, matrix3, sizeof(matrix3) );
	LTR_VERIFY( ltr_MeshAddInstance( mesh1, &m1inst3 ) );
	
	// LIGHTS
	ltr_LightInfo lights[] =
	{
		{ LTR_LT_POINT, { -0.4f, -0.4f, 1 }, {0,0,0}, {0,0,0}, { 0.9f, 0.1f, 0.05f }, 4.0f, 1.0f, 0.1f, 9 },
	};
	for( size_t lt = 0; lt < sizeof(lights)/sizeof(lights[0]); ++lt )
		ltr_LightAdd( scene, &lights[ lt ] );
	
	// --- DO WORK ---
	DOWORK( scene );
	
	// --- RETURN OUTPUT ---
	ltr_WorkOutput wout;
	ltr_WorkOutputInfo woutinfo;
	ltr_GetWorkOutputInfo( scene, &woutinfo );
	for( u32 lm = 0; lm < woutinfo.lightmap_count; ++lm )
	{
		LTR_VERIFY( ltr_GetWorkOutput( scene, lm, &wout ) );
		char namebfr[ 64 ];
		sprintf( namebfr, "%d.tga", (int) wout.uid );
		dumpimg( namebfr, wout.lightmap_rgb, wout.width, wout.height );
	}
	
	ltr_DestroyScene( scene );
}

void testfunc_mesh1()
{
	puts( "LIGHTER test [mesh1]" );
	
	ltr_Scene* scene = ltr_CreateScene();
	
	ltr_Config cfg;
	ltr_GetConfig( &cfg, scene );
	cfg.ao_distance = 2;
	ltr_SetConfig( scene, &cfg );
	
	// MESH 1
	ltr_Mesh* mesh1 = ltr_CreateMesh( scene, "mesh1", strlen("mesh1") );
	testmesh M;
	LTR_VERIFY( (int) M.Load( "bin/test-mesh.data" ) );
	
	// MESH 1 PART 1
	ltr_MeshPartInfo m1part1 =
	{
		&M.positions[0], &M.normals[0], &M.texcoords[0], &M.texcoords[0],
		sizeof(float)*3, sizeof(float)*3, sizeof(float)*2, sizeof(float)*2,
		&M.indices[0], (u32) M.positions.size(), (u32) M.indices.size(), 1
	};
	LTR_VERIFY( ltr_MeshAddPart( mesh1, &m1part1 ) );
	
	// MESH 1 INSTANCE 1
	ltr_MeshInstanceInfo m1inst1;
	memset( &m1inst1, 0, sizeof(m1inst1) );
	m1inst1.importance = 1;
	m1inst1.shadow = 1;
	float matrix1[16] = { 1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
	memcpy( m1inst1.matrix, matrix1, sizeof(matrix1) );
	LTR_VERIFY( ltr_MeshAddInstance( mesh1, &m1inst1 ) );
	
	// LIGHTS
	ltr_LightInfo lights[] =
	{
		{ LTR_LT_POINT, { -2.18f, -4.04f, 1.40f }, {0,0,0}, {0,0,0}, { 0.9f, 0.7f, 0.5f }, 16.0f, 1.0f, 0.1f, 5 },
		{ LTR_LT_POINT, { 2.18f, 4.04f, 1.40f }, {0,0,0}, {0,0,0}, { 0.5f, 0.7f, 0.9f }, 16.0f, 1.0f, 0.1f, 5 },
		{ LTR_LT_SPOT, { 0, 0, 1.60f }, {0,0,-1}, {1,0,0}, { 0.7f, 0.1f, 0.05f }, 16.0f, 1.0f, 0.1f, 5, 45.0f, 25.0f, 0.5f },
	};
	for( size_t lt = 0; lt < sizeof(lights)/sizeof(lights[0]); ++lt )
		ltr_LightAdd( scene, &lights[ lt ] );
	
	// --- DO WORK ---
	DOWORK( scene );
	
	// --- RETURN OUTPUT ---
	ltr_WorkOutput wout;
	ltr_WorkOutputInfo woutinfo;
	ltr_GetWorkOutputInfo( scene, &woutinfo );
	for( u32 lm = 0; lm < woutinfo.lightmap_count; ++lm )
	{
		LTR_VERIFY( ltr_GetWorkOutput( scene, lm, &wout ) );
		char namebfr[ 64 ];
		sprintf( namebfr, "%d.tga", (int) wout.uid );
		dumpimg( namebfr, wout.lightmap_rgb, wout.width, wout.height );
	}
	
	ltr_DestroyScene( scene );
}

void testfunc_mesh2()
{
	puts( "LIGHTER test [mesh2]" );
	
	ltr_Scene* scene = ltr_CreateScene();
	
	ltr_Config cfg;
	ltr_GetConfig( &cfg, scene );
	cfg.ao_distance = 2;
	cfg.global_size_factor = 8;
	cfg.blur_size = 1;
//	cfg.ds2x = 1;
	ltr_SetConfig( scene, &cfg );
	
	// MESH 1
	ltr_Mesh* mesh1 = ltr_CreateMesh( scene, "mesh1", strlen("mesh1") );
	testmesh M;
	LTR_VERIFY( (int) M.Load( "bin/test-set2-mesh1.data" ) );
	
	// MESH 1 PART 1
	ltr_MeshPartInfo m1part1 =
	{
		&M.positions[0], &M.normals[0], &M.texcoords[0], &M.texcoords[0],
		sizeof(float)*3, sizeof(float)*3, sizeof(float)*2, sizeof(float)*2,
		&M.indices[0], (u32) M.positions.size(), (u32) M.indices.size(), 1
	};
	LTR_VERIFY( ltr_MeshAddPart( mesh1, &m1part1 ) );
	
	// MESH 1 INSTANCE 1
	ltr_MeshInstanceInfo m1inst1;
	memset( &m1inst1, 0, sizeof(m1inst1) );
	m1inst1.importance = 1;
	m1inst1.shadow = 1;
	float matrix1[16] = { 1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
	memcpy( m1inst1.matrix, matrix1, sizeof(matrix1) );
	LTR_VERIFY( ltr_MeshAddInstance( mesh1, &m1inst1 ) );
	
	// MESH 2
	ltr_Mesh* mesh2 = ltr_CreateMesh( scene, "mesh2", strlen("mesh2") );
	testmesh M2;
	LTR_VERIFY( (int) M2.Load( "bin/test-set2-mesh2.data" ) );
	
	// MESH 2 PART 1
	ltr_MeshPartInfo m2part1 =
	{
		&M2.positions[0], &M2.normals[0], &M2.texcoords[0], &M2.texcoords[0],
		sizeof(float)*3, sizeof(float)*3, sizeof(float)*2, sizeof(float)*2,
		&M2.indices[0], (u32) M2.positions.size(), (u32) M2.indices.size(), 1
	};
	LTR_VERIFY( ltr_MeshAddPart( mesh2, &m2part1 ) );
	
	// MESH 2 INSTANCE 1
	ltr_MeshInstanceInfo m2inst1;
	memset( &m2inst1, 0, sizeof(m2inst1) );
	m2inst1.importance = 0.3f;
	m2inst1.shadow = 1;
	float matrix2[16] = { 1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
	memcpy( m2inst1.matrix, matrix2, sizeof(matrix2) );
	LTR_VERIFY( ltr_MeshAddInstance( mesh2, &m2inst1 ) );
	
	// LIGHTS
	ltr_LightInfo lights[] =
	{
		{ LTR_LT_POINT, { -2.18f, -4.04f, 1.40f }, {0,0,0}, {0,0,0}, { 0.9f, 0.7f, 0.5f }, 16.0f, 1.0f, 0.2f, 5,  0, 0, 0 },
		{ LTR_LT_POINT, { 2.18f, 4.04f, 1.40f }, {0,0,0}, {0,0,0}, { 0.5f, 0.7f, 0.9f }, 16.0f, 1.0f, 0.2f, 5,  0, 0, 0 },
		{ LTR_LT_SPOT, { 0, 0, 1.60f }, {0,0,-1}, {1,0,0}, { 0.7f, 0.1f, 0.05f }, 16.0f, 1.0f, 0.2f, 5, 45.0f, 25.0f, 0.5f },
		{ LTR_LT_DIRECT, {0,0,0}, { 1, 1, 1 }, {0,0,0}, { 0.6f, 0.55f, 0.5f }, 1000.0f, 1.0f, 0.2f, 1,  0, 0, 0 },
	};
	for( size_t lt = 0; lt < sizeof(lights)/sizeof(lights[0]); ++lt )
		ltr_LightAdd( scene, &lights[ lt ] );
	
	// --- DO WORK ---
	DOWORK( scene );
	
	// --- RETURN OUTPUT ---
	ltr_WorkOutput wout;
	ltr_WorkOutputInfo woutinfo;
	ltr_GetWorkOutputInfo( scene, &woutinfo );
	for( u32 lm = 0; lm < woutinfo.lightmap_count; ++lm )
	{
		LTR_VERIFY( ltr_GetWorkOutput( scene, lm, &wout ) );
		char namebfr[ 64 ];
		sprintf( namebfr, "%d.tga", (int) wout.uid );
		dumpimg( namebfr, wout.lightmap_rgb, wout.width, wout.height );
	}
	
	ltr_DestroyScene( scene );
}

LTRBOOL rad1_samplefunc( ltr_Config* cfg, ltr_SampleRequest* req )
{
//	printf( "T1[ %g %g ] T2[ %g %g ]\n", req->tex0u, req->tex0v, req->tex1u, req->tex1v );
	if( req->position[0] <= -3 && Vec3Dot( Vec3::CreateFromPtr( req->normal ), Vec3::Create( 1, 0, 0 ) ) > 0.1f )
	{
		req->out_diffuse_color[0] = 0.5f;
		req->out_diffuse_color[1] = 0.05f;
		req->out_diffuse_color[2] = 0.02f;
	}
	return 1;
}

void testfunc_rad1()
{
	puts( "LIGHTER test [rad1]" );
	
	ltr_Scene* scene = ltr_CreateScene();
	
	ltr_Config cfg;
	ltr_GetConfig( &cfg, scene );
	cfg.ao_distance = 2;
	cfg.bounce_count = 2;
	cfg.sample_fn = rad1_samplefunc;
	ltr_SetConfig( scene, &cfg );
	
	// MESH 1
	ltr_Mesh* mesh1 = ltr_CreateMesh( scene, "mesh1", strlen("mesh1") );
	testmesh M;
	LTR_VERIFY( (int) M.Load( "bin/test-set2-mesh1.data" ) );
	
	// MESH 1 PART 1
	ltr_MeshPartInfo m1part1 =
	{
		&M.positions[0], &M.normals[0], &M.texcoords[0], &M.texcoords[0],
		sizeof(float)*3, sizeof(float)*3, sizeof(float)*2, sizeof(float)*2,
		&M.indices[0], (u32) M.positions.size(), (u32) M.indices.size(), 1
	};
	LTR_VERIFY( ltr_MeshAddPart( mesh1, &m1part1 ) );
	
	// MESH 1 INSTANCE 1
	ltr_MeshInstanceInfo m1inst1;
	memset( &m1inst1, 0, sizeof(m1inst1) );
	m1inst1.importance = 0.4f;
	m1inst1.shadow = 1;
	float matrix1[16] = { 1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
	memcpy( m1inst1.matrix, matrix1, sizeof(matrix1) );
	LTR_VERIFY( ltr_MeshAddInstance( mesh1, &m1inst1 ) );
	
	// LIGHTS
	ltr_LightInfo lights[] =
	{
		{ LTR_LT_POINT, { 10.18f, 0, 1.40f }, {0,0,0}, {0,0,0}, { 0.99f, 0.97f, 0.95f }, 24.0f, 1.0f, 0.1f, 5 },
//		{ LTR_LT_POINT, { 2.18f, 4.04f, 1.40f }, {0,0,0}, {0,0,0}, { 0.5f, 0.7f, 0.9f }, 16.0f, 1.0f, 0.1f, 5 },
//		{ LTR_LT_SPOT, { 0, 0, 1.60f }, {0,0,-1}, {1,0,0}, { 0.7f, 0.1f, 0.05f }, 16.0f, 1.0f, 0.1f, 5, 45.0f, 25.0f, 0.5f },
	};
	for( size_t lt = 0; lt < sizeof(lights)/sizeof(lights[0]); ++lt )
		ltr_LightAdd( scene, &lights[ lt ] );
	
	// --- DO WORK ---
	DOWORK( scene );
	
	// --- RETURN OUTPUT ---
	ltr_WorkOutput wout;
	ltr_WorkOutputInfo woutinfo;
	ltr_GetWorkOutputInfo( scene, &woutinfo );
	for( u32 lm = 0; lm < woutinfo.lightmap_count; ++lm )
	{
		LTR_VERIFY( ltr_GetWorkOutput( scene, lm, &wout ) );
		char namebfr[ 64 ];
		sprintf( namebfr, "%d.tga", (int) wout.uid );
		dumpimg( namebfr, wout.lightmap_rgb, wout.width, wout.height );
	}
	
	ltr_DestroyScene( scene );
}


const char* testname = "mesh2";

int main( int argc, char* argv[] )
{
	if( argc > 1 )
		testname = argv[1];
	
	if( strcmp( testname, "basic" ) == 0 )
		testfunc_basic();
	else if( strcmp( testname, "mesh1" ) == 0 )
		testfunc_mesh1();
	else if( strcmp( testname, "mesh2" ) == 0 )
		testfunc_mesh2();
	else if( strcmp( testname, "rad1" ) == 0 )
		testfunc_rad1();
	else
		printf( "TEST NOT FOUND: %s\n", testname );
	
	return 0;
}
