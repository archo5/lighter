
// LIGHtmap
// TExture
// Renderer

#pragma once

#ifdef __cplusplus
extern "C" {
#endif


#if defined(_WIN32) && !defined(LTR_STATIC)
#ifdef LTRBUILD
#define LTRAPI __declspec(dllexport)
#else
#define LTRAPI __declspec(dllimport)
#endif
#else
#define LTRAPI
#endif


#define LTRBOOL int
#define LTRCODE int

#include <stdlib.h>
#ifdef _MSC_VER
#  pragma warning( disable: 4996 )
#  if _MSC_VER >= 1600
#    include <stdint.h>
#  else
typedef __int8 int8_t;
typedef __int16 int16_t;
typedef __int32 int32_t;
typedef __int64 int64_t;
typedef unsigned __int8 uint8_t;
typedef unsigned __int16 uint16_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;
#  endif
#else
#  include <inttypes.h>
#endif
typedef int32_t i32;
typedef uint32_t u32;

typedef float ltr_VEC3[3];
typedef float ltr_VEC4[4];
typedef ltr_VEC4 ltr_MAT4[4];


#define LTR_VEC3_SET( v, x, y, z ) do{ (v)[0]=(x); (v)[1]=(y); (v)[2]=(z); }while(0)


#define LTRC_SUCCESS 0
#define LTRC_ATEND   1

#define LTR_LT_POINT  1
#define LTR_LT_SPOT   2
#define LTR_LT_DIRECT 3


typedef struct ltr_MeshPartInfo
{
	const void* positions_f3;
	const void* normals_f3;
	const void* texcoords1_f2;
	const void* texcoords2_f2;
	u32 stride_positions;
	u32 stride_normals;
	u32 stride_texcoords1;
	u32 stride_texcoords2;
	const u32* indices;
	u32 vertex_count;
	u32 index_count;
	int shadow;
}
ltr_MeshPartInfo;

typedef struct ltr_MeshInstanceInfo
{
	ltr_MAT4 matrix;
	float importance;
	int shadow;
	const char* ident;
	size_t ident_size;
}
ltr_MeshInstanceInfo;

typedef struct ltr_LightInfo
{
	u32 type;
	ltr_VEC3 position;
	ltr_VEC3 direction;
	ltr_VEC3 up_direction;
	ltr_VEC3 color_rgb;
	float range;
	float power;
	float light_radius;
	int shadow_sample_count;
	float spot_angle_out;
	float spot_angle_in;
	float spot_curve;
}
ltr_LightInfo;

typedef struct ltr_SampleInfo
{
	u32 id;
	ltr_VEC3 position;
	ltr_VEC3 normal;
	ltr_VEC3 out_color;
}
ltr_SampleInfo;

typedef struct ltr_SampleRequest
{
	ltr_VEC3 position;
	ltr_VEC3 normal;
	float tex0u, tex0v;
	float tex1u, tex1v;
	uint32_t part_id;
	
	const char* mesh_ident;
	size_t mesh_ident_size;
	const char* inst_ident;
	size_t inst_ident_size;
	
	ltr_VEC3 out_diffuse_color;
	ltr_VEC3 out_emissive_color;
}
ltr_SampleRequest;

typedef struct ltr_Config
{
	// callbacks
	void* userdata;
	LTRBOOL (*size_fn)
	(
		ltr_Config* /* config */,
		const char* /* mesh_ident */,
		size_t /* mesh_ident_size */,
		const char* /* inst_ident */,
		size_t /* inst_ident_size */,
		float /* computed_surface_area */,
		float /* inst_importance */,
		u32[2] /* out_width, out_height */
	);
	// limits & factors
	int max_num_threads;
	size_t max_tree_memory;
	u32 max_lightmap_size;
	u32 default_width;
	u32 default_height;
	float global_size_factor;
	float max_correct_dist;
	float max_correct_angle;
	// lightmap generation
	ltr_VEC3 clear_color; /* the color used outside triangles */
	ltr_VEC3 ambient_color; /* the base color used inside triangles */
	// RADIOSITY:
	int bounce_count;
	LTRBOOL (*sample_fn)
	(
		ltr_Config* /* config */,
		ltr_SampleRequest* /* req */
	);
	// AMBIENT OCCLUSION effect:
	// ao_factor = 1 - max( raytrace_distance / ao_distance, 1 )
	// - AO FACTORS are accumulated -
	// ao_factor = pow( min( ao_factor * ao_multiplier, 1 ), ao_falloff )
	// IF ao_effect >= 0: color = color_orig * ( 1 - ao_effect * ao_factor ) + color_ao * ao_factor;
	// IF ao_effect <  0: color = lerp( color_orig, lerp( color_ao, color_orig * color_ao, -ao_effect ), ao_factor )
	float ao_distance;
	float ao_multiplier;
	float ao_falloff;
	float ao_effect;
	float ao_divergence; /* hemisphere [0] -to- full sphere [1] direction picking */
	ltr_VEC3 ao_color_rgb;
	int ao_num_samples;
	// GAUSSIAN BLUR effect:
	float blur_size;
	// DOWNSAMPLE 2X effect:
	int ds2x;
	// NORMAL/FOCUS data generation:
	int generate_normalmap_data;
}
ltr_Config;

LTRAPI LTRBOOL ltr_DefaultSizeFunc
(
	ltr_Config* config,
	const char* mesh_ident,
	size_t mesh_ident_size,
	const char* inst_ident,
	size_t inst_ident_size,
	float computed_surface_area,
	float inst_importance,
	u32 out_size[2]
);

typedef struct ltr_WorkOutputInfo
{
	u32 lightmap_count;
	u32 sample_count;
	ltr_SampleInfo* samples;
}
ltr_WorkOutputInfo;

typedef struct ltr_WorkOutput
{
	u32 uid;
	const char* mesh_ident;
	size_t mesh_ident_size;
	const char* inst_ident;
	size_t inst_ident_size;
	float* lightmap_rgb;
	float* normals_xyzf;
	u32 width;
	u32 height;
}
ltr_WorkOutput;

typedef struct ltr_Mesh ltr_Mesh;
typedef struct ltr_Scene ltr_Scene;

typedef struct ltr_WorkStatus
{
	float completion;
	const char* stage;
}
ltr_WorkStatus;


// - manage scene
LTRAPI ltr_Scene* ltr_CreateScene();
LTRAPI void ltr_DestroyScene( ltr_Scene* scene );
LTRAPI void ltr_Start( ltr_Scene* scene );
LTRAPI void ltr_Abort( ltr_Scene* scene );
LTRAPI LTRBOOL ltr_GetStatus( ltr_Scene* scene, ltr_WorkStatus* wsout );
LTRAPI void ltr_Sleep( int ms );
LTRAPI void ltr_GetConfig( ltr_Config* cfg, ltr_Scene* opt_scene );
LTRAPI LTRCODE ltr_SetConfig( ltr_Scene* scene, ltr_Config* cfg );

// - data input
LTRAPI ltr_Mesh* ltr_CreateMesh( ltr_Scene* scene, const char* ident, size_t ident_size );
LTRAPI LTRBOOL ltr_MeshAddPart( ltr_Mesh* mesh, ltr_MeshPartInfo* mpinfo );
LTRAPI LTRBOOL ltr_MeshAddInstance( ltr_Mesh* mesh, ltr_MeshInstanceInfo* mii );
LTRAPI void ltr_LightAdd( ltr_Scene* scene, ltr_LightInfo* li );
LTRAPI void ltr_SampleAdd( ltr_Scene* scene, ltr_SampleInfo* si );

// - data output
LTRAPI void ltr_GetWorkOutputInfo( ltr_Scene* scene, ltr_WorkOutputInfo* woutinfo );
LTRAPI LTRBOOL ltr_GetWorkOutput( ltr_Scene* scene, u32 which, ltr_WorkOutput* wout );

// - misc.
LTRAPI u32 ltr_NextPowerOfTwo( u32 x );


#ifdef __cplusplus
}
#endif

