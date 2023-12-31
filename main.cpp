#include <Novice.h>
#include"MathFunction.h"
#include<stdint.h>

const char kWindowTitle[] = "LE2B_28_ヨシダアイリ_クロス実績と3d描画";

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	const int kWindowWidth = 1280;
	const int kWindowHeight = 720;
	Novice::Initialize(kWindowTitle, kWindowWidth, kWindowHeight);

	// キー入力結果を受け取る箱
	char keys[256] = {0};
	char preKeys[256] = {0};

	//クロス積用
	Vector3 v1{ 1.2f,-3.9f,2.5f };
	Vector3 v2{ 2.8f,0.4f,-1.3f };

	//三角形描画
	Vector3 rotate{};
	Vector3 translate{};
	Vector3 cameraPosition{ 0.0f,0.0f,-10.0f };
	Vector3 kLocalVertics[3] =
	{
		/*{0.0f,0.3f,1.0f},
		{0.6f,-0.6f,1.0f},
		{-0.6f,-0.6f,1.0f},*/

		{0.0f,1.0f,0.0f},
		{-1.0f,0.0f,0.0f},
		{1.0f,0.0f,0.0f}
	};

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		//クロス積
		Vector3 cross = Cross(v1, v2);

		//前後移動
		if (keys[DIK_W])
		{
			translate.z -= 0.1f;
		}
		if (keys[DIK_S])
		{
			translate.z += 0.1f;
		}
		//左右移動
		if (keys[DIK_A])
		{
			translate.x -= 0.1f;
		}
		if (keys[DIK_D])
		{
			translate.x += 0.1f;
		}

		//回転

		rotate.y += 0.1f;

		//各種行列の計算
		Matrix4x4 worldMatrix = MakeAffinMatrix({ 1.0f,1.0f,1.0f }, rotate, translate);
		Matrix4x4 cameraMatrix = MakeAffinMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, cameraPosition);
		Matrix4x4 viewMatrix = Inverse(cameraMatrix);
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
		Matrix4x4 worldViewProjectionmatirx = Multiply(worldMatrix, Multiply(viewMatrix, projectionMatrix));
		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);
		Vector3 screenVertices[3];
		for (uint32_t i = 0; i < 3; ++i)
		{
			Vector3 ndcVertex = Transform(kLocalVertics[i], worldViewProjectionmatirx);
			screenVertices[i] = Transform(ndcVertex, viewportMatrix);
		}

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		//クロス積
		VectorScreenPrintf(0, 0, cross, "Cross");

		//三角形の描画
		Novice::DrawTriangle
		(
			int(screenVertices[0].x), int(screenVertices[0].y),
			int(screenVertices[1].x), int(screenVertices[1].y),
			int(screenVertices[2].x), int(screenVertices[2].y),
			RED, kFillModeSolid
		);

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
