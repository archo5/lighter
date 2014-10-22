
// LIGHtmap
// TExture
// Renderer

#pragma once

#ifdef __cplusplus
extern "C" {
#endif


#ifdef LTRBUILD
#define LTRAPI __declspec(dllexport)
#else
#define LTRAPI __declspec(dllimport)
#endif


#define LTRBOOL int
#define LTRCODE int

#include <stdlib.h>
#include <inttypes.h>
typedef int32_t i32;
typedef uint32_t u32;

typedef float ltr_VEC3[3];
typedef float ltr_VEC4[4];
typedef ltr_VEC4 ltr_MAT4[4];


#define LTRC_SUCCESS 0
#define LTRC_ATEND   1

#define LTR_LT_POINT  1
#define LTR_LT_SPOT   2
#define LTR_LT_DIRECT 3


typedef struct ltr_WorkInfo
{
	i32 type;
	i32 part;
	const char* stage;
	i32 item_count;
}
ltr_WorkInfo;

typedef struct ltr_MeshPartInfo
{
	void* positions_f3;
	void* normals_f3;
	void* texcoords1_f2;
	void* texcoords2_f2;
	u32 stride_positions;
	u32 stride_normals;
	u32 stride_texcoords1;
	u32 stride_texcoords2;
	u32* indices;
	u32 vertex_count;
	u32 index_count;
	int tristrip;
}
ltr_MeshPartInfo;

typedef struct ltr_MeshInstanceInfo
{
	ltr_MAT4 matrix;
	float importance;
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
}
ltr_LightInfo;

typedef struct ltr_WorkOutputInfo
{
	i32 lightmap_count;
}
ltr_WorkOutputInfo;

typedef struct ltr_WorkOutput
{
	u32 uid;
	const char* ident;
	size_t ident_size;
	float* lightmap_rgb;
	u32 width;
	u32 height;
}
ltr_WorkOutput;

typedef struct ltr_Mesh ltr_Mesh;
typedef struct ltr_Scene ltr_Scene;


// - manage scene
LTRAPI ltr_Scene* ltr_CreateScene();
LTRAPI void ltr_DestroyScene( ltr_Scene* scene );
LTRAPI LTRCODE ltr_DoWork( ltr_Scene* scene, ltr_WorkInfo* info );

// - data input
LTRAPI ltr_Mesh* ltr_CreateMesh( ltr_Scene* scene, const char* ident, size_t ident_size );
LTRAPI LTRBOOL ltr_MeshAddPart( ltr_Mesh* mesh, ltr_MeshPartInfo* mpinfo );
LTRAPI LTRBOOL ltr_MeshAddInstance( ltr_Mesh* mesh, ltr_MeshInstanceInfo* mii );
LTRAPI void ltr_LightAdd( ltr_Scene* scene, ltr_LightInfo* li );

// - data output
LTRAPI void ltr_GetWorkOutputInfo( ltr_Scene* scene, ltr_WorkOutputInfo* woutinfo );
LTRAPI LTRBOOL ltr_GetWorkOutput( ltr_Scene* scene, i32 which, ltr_WorkOutput* wout );


#ifdef __cplusplus
}
#endif

