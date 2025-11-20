#include <math.h>

#if defined __linux__
#include "/home/codeleaded/System/Static/Library/ConsoleEngine.h"
#include "/home/codeleaded/System/Static/Container/Vector.h"
#elif defined _WINE
#include "/home/codeleaded/System/Static/Library/ConsoleEngine.h"
#include "/home/codeleaded/System/Static/Container/Vector.h"
#elif defined _WIN32
#include "F:/home/codeleaded/System/Static/Library/ConsoleEngine.h"
#include "F:/home/codeleaded/System/Static/Container/Vector.h"
#elif defined(__APPLE__)
#error "Apple not supported!"
#else
#error "Platform not supported!"
#endif

typedef struct Vec3D{
	float x, y, z;
} Vec3D;

typedef struct Tri3D{
	Vec3D p[3];
} Tri3D;

typedef struct mesh{
	Vector tris;
} mesh;

typedef struct M4x4D{
	float m[4][4];
} M4x4D;

mesh meshCube;
M4x4D matProj;

float fTheta;

void MultiplyMatrixVector(Vec3D i, Vec3D* o, M4x4D m){
	o->x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
	o->y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
	o->z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];
	float w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];
	if (w != 0.0f)
	{
		o->x /= w; o->y /= w; o->z /= w;
	}
}

void C_Setup(Console* c){
    meshCube.tris = Vector_New(sizeof(Tri3D));
    
    // SOUTH                                                     
    Tri3D tri = {0.0f, 0.0f, 0.0f,    0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 0.0f};
    Vector_Push(&meshCube.tris,&tri);
    tri = (Tri3D){0.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 0.0f, 0.0f};
    Vector_Push(&meshCube.tris,&tri);
    
    //EAST                                                     
    tri = (Tri3D){ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f };
    Vector_Push(&meshCube.tris,&tri);
    tri = (Tri3D){ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 0.0f, 1.0f };
    Vector_Push(&meshCube.tris,&tri);
    
    // NORTH                                                     
    tri = (Tri3D){ 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f };
    Vector_Push(&meshCube.tris,&tri);
    tri = (Tri3D){ 1.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f };
    Vector_Push(&meshCube.tris,&tri);
    
    // WEST                                                      
    tri = (Tri3D){ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f };
    Vector_Push(&meshCube.tris,&tri);
    tri = (Tri3D){ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 0.0f,    0.0f, 0.0f, 0.0f };
    Vector_Push(&meshCube.tris,&tri);
    
    // TOP                                                       
    tri = (Tri3D){ 0.0f, 1.0f, 0.0f,    0.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f };
    Vector_Push(&meshCube.tris,&tri);
    tri = (Tri3D){ 0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 0.0f };
    Vector_Push(&meshCube.tris,&tri);
    
    // BOTTOM                                                    
    tri = (Tri3D){ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f };
    Vector_Push(&meshCube.tris,&tri);
    tri = (Tri3D){ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f,    1.0f, 0.0f, 0.0f };
    Vector_Push(&meshCube.tris,&tri);


	float fNear = 0.1f;
	float fFar = 1000.0f;
	float fFov = 90.0f;
	float fAspectRatio = (float)Console_Height(c) / (float)Console_Width(c);
	float fFovRad = 1.0f / tanf(fFov * 0.5f / 180.0f * 3.14159f);
	memset(matProj.m[0],0,4 * sizeof(float));
    memset(matProj.m[1],0,4 * sizeof(float));
    memset(matProj.m[2],0,4 * sizeof(float));
    memset(matProj.m[3],0,4 * sizeof(float));
    matProj.m[0][0] = fAspectRatio * fFovRad;
	matProj.m[1][1] = fFovRad;
	matProj.m[2][2] = fFar / (fFar - fNear);
	matProj.m[3][2] = (-fFar * fNear) / (fFar - fNear);
	matProj.m[2][3] = 1.0f;
	matProj.m[3][3] = 0.0f;
}
void C_Update(Console* c){
    // Clear Screen
	Console_Clear(c,CPIXEL_SOLID,FG_BLACK);
	
	// Set up rotation matrices
	M4x4D matRotZ, matRotX;
	fTheta += 1.0f * c->ElapsedTime;
	
    // Rotation Z
    memset(matRotZ.m[0],0,4 * sizeof(float));
    memset(matRotZ.m[1],0,4 * sizeof(float));
    memset(matRotZ.m[2],0,4 * sizeof(float));
    memset(matRotZ.m[3],0,4 * sizeof(float));
	matRotZ.m[0][0] = cosf(fTheta);
	matRotZ.m[0][1] = sinf(fTheta);
	matRotZ.m[1][0] = -sinf(fTheta);
	matRotZ.m[1][1] = cosf(fTheta);
	matRotZ.m[2][2] = 1;
	matRotZ.m[3][3] = 1;

	// Rotation X
    memset(matRotX.m[0],0,4 * sizeof(float));
    memset(matRotX.m[1],0,4 * sizeof(float));
    memset(matRotX.m[2],0,4 * sizeof(float));
    memset(matRotX.m[3],0,4 * sizeof(float));
	matRotX.m[0][0] = 1;
	matRotX.m[1][1] = cosf(fTheta * 0.5f);
	matRotX.m[1][2] = sinf(fTheta * 0.5f);
	matRotX.m[2][1] = -sinf(fTheta * 0.5f);
	matRotX.m[2][2] = cosf(fTheta * 0.5f);
	matRotX.m[3][3] = 1;
	
	// Draw Triangles
	for (int i = 0;i<meshCube.tris.size;i++){
        Tri3D tri = *(Tri3D*)Vector_Get(&meshCube.tris,i);
		Tri3D triProjected, triTranslated, triRotatedZ, triRotatedZX;
		// Rotate in Z-Axis
		MultiplyMatrixVector(tri.p[0], &triRotatedZ.p[0], matRotZ);
		MultiplyMatrixVector(tri.p[1], &triRotatedZ.p[1], matRotZ);
		MultiplyMatrixVector(tri.p[2], &triRotatedZ.p[2], matRotZ);
		// Rotate in X-Axis
		MultiplyMatrixVector(triRotatedZ.p[0], &triRotatedZX.p[0], matRotX);
		MultiplyMatrixVector(triRotatedZ.p[1], &triRotatedZX.p[1], matRotX);
		MultiplyMatrixVector(triRotatedZ.p[2], &triRotatedZX.p[2], matRotX);
		// Offset into the screen
		triTranslated = triRotatedZX;
		triTranslated.p[0].z = triRotatedZX.p[0].z + 3.0f;
		triTranslated.p[1].z = triRotatedZX.p[1].z + 3.0f;
		triTranslated.p[2].z = triRotatedZX.p[2].z + 3.0f;
		// Project Tri3Ds from 3D --> 2D
		MultiplyMatrixVector(triTranslated.p[0], &triProjected.p[0], matProj);
		MultiplyMatrixVector(triTranslated.p[1], &triProjected.p[1], matProj);
		MultiplyMatrixVector(triTranslated.p[2], &triProjected.p[2], matProj);
		// Scale into view
		triProjected.p[0].x += 1.0f; triProjected.p[0].y += 1.0f;
		triProjected.p[1].x += 1.0f; triProjected.p[1].y += 1.0f;
		triProjected.p[2].x += 1.0f; triProjected.p[2].y += 1.0f;
		triProjected.p[0].x *= 0.5f * (float)Console_Width(c);
		triProjected.p[0].y *= 0.5f * (float)Console_Height(c);
		triProjected.p[1].x *= 0.5f * (float)Console_Width(c);
		triProjected.p[1].y *= 0.5f * (float)Console_Height(c);
		triProjected.p[2].x *= 0.5f * (float)Console_Width(c);
		triProjected.p[2].y *= 0.5f * (float)Console_Height(c);
		// Rasterize Tri3D
		Console_RenderTriangleWire(
			c,
			(Vec2){triProjected.p[0].x,triProjected.p[0].y},
			(Vec2){triProjected.p[1].x,triProjected.p[1].y},
			(Vec2){triProjected.p[2].x,triProjected.p[2].y},
			(CPixel){CPIXEL_SOLID, FG_WHITE},
			1.0f
		);
	}
}
void C_Delete(Console* c){
    Vector_Free(&meshCube.tris);
}

int main(){
    Console c;
    if(Console_Create(&c,L"3D Engine",200,150,8,8,C_Setup,C_Update,C_Delete)){
        Console_Start(&c);
    }
    return 0;
}