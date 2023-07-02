﻿#pragma once
#include"Vector3.h"
#include"Matrix4x4.h"

//1.行列の加法
Matrix4x4 Add(const Matrix4x4& m1, const Matrix4x4& m2);
//2.行列の減法
Matrix4x4 Subtract(const Matrix4x4& m1, const Matrix4x4& m2);
//3.行列の積
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2);
//4.逆行列
Matrix4x4 Inverse(const Matrix4x4& m);
////5.転置行列
Matrix4x4 Transpose(const Matrix4x4& m);
////6.単位行列の作成
Matrix4x4 MakeIdentity4x4();

//1.平行移動行列
Matrix4x4 MakeTranselateMatrix(const Vector3& transelate);
//2.拡大縮小行列
Matrix4x4 MakeScaleMatrix(const Vector3& scale);
//3.座標変換
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix);

//1.x軸回転行列
Matrix4x4 MakeRotateXMatrix(float radian);
//2.y軸回転行列
Matrix4x4 MakeRotateYMatrix(float radian);
//3.z軸回転行列
Matrix4x4 MakeRotateZMatrix(float radian);

//3次元アフィン変換行列
Matrix4x4 MakeAffinMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate);

//1.透視投影行列
Matrix4x4 MakePerspectiveFovMatrix(float forY, float aspectRatio, float nearClip, float farClip);
//2.正射影行列
Matrix4x4 MakeOrthographicmatrix(float left, float top, float right, float bottom, float nearClip, float farClip);

//3.ビューポート変換行列
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth);

//クロス積
Vector3 Cross(const Vector3& v1, const Vector3& v2);

//描画処理
//四次元
static const int kRowHeight = 20;
static const int kColumnWidth = 60;
void MatrixScreenPrintf(int x, int y, const Matrix4x4& matrix, const char* label);
//三次元
void VectorScreenPrintf(int x, int y, const Vector3& vector, const char* label);