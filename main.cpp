
#include <Novice.h>
#include <cmath>
#include <cassert>
#include <imgui.h>
#include <algorithm>

const char kWindowTitle[] = "LE2B_26";

typedef struct Vector3 {
	float x;
	float y;
	float z;
}Vector3;

typedef struct Matrix4x4 {
	float m[4][4];
}Matrix4x4;

struct Line {
	Vector3 origin;
	Vector3 diff;
};

struct Ray {
	Vector3 origin;
	Vector3 diff;
};

struct Segment {
	Vector3 origin;
	Vector3 diff;
};

struct Sphere {
	Vector3 center;
	float radius;
};

struct Plane {
	Vector3 normal;
	float distance;
};

struct Triangle {
	Vector3 vertices[3];
};

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	Segment segment{ {-2.0f, -1.0f, 0.0f}, {3.0f, 2.0f, 2.0f} };
	Line line{ {-2.0f, -1.0f, 0.0f}, {3.0f, 2.0f, 2.0f} };
	Plane plane{ {0.0f,1.0f,0.0f}, 0.0f };
	Triangle triangle{ {{-1.0f,0.0f, 0.0f}, {1.0f,0.0f,0.0f}, {0.0f, 2.0f, 0.0f}} };
	Vector3 point{ -1.5f, 0.6f, 0.6f };
	Vector3 cameraTranslate{ -1.0f, 1.9f,-6.49f };
	Vector3 cameraRotate{ 0.0f, 3.0f,-15.7f };
	Vector3 rotate{ 0.0f, 0.0f, 0.0f };
	Vector3 translate{ 0.0f, 0.0f, 5.57f };

	Vector3 Transform(const Vector3 & vector, const Matrix4x4 & matrix);
	Vector3 Add(const Vector3 & v1, const Vector3 & v2);
	Vector3 Subtract(const Vector3 & v1, const Vector3 & v2);
	Matrix4x4 Multiply(const Matrix4x4 & m1, const Matrix4x4 & m2);
	Matrix4x4 Inverse(const Matrix4x4 & m1);
	Vector3 Normalize(Vector3 v1);
	Vector3 Project(const Vector3 & v1, const Vector3 & v2);
	Vector3 ClosestPoint(const Vector3 & point, const Segment & segment);
	Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth);
	Matrix4x4 MakeAffineMatrix(const Vector3 & scale, const Vector3 & rotate, const Vector3 & translate);
	Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip);
	bool IsCollision(const Sphere & s1, const Sphere & s2);
	bool IsCollision(const Sphere & s, const Plane & p);
	bool IsCollision(const Segment & segment, const Plane & plane);
	bool IsCollision(const Triangle & triangle, const Segment & segment);
	void DrawSphere(const Sphere & sphere, const Matrix4x4 & viewProjectionMatrix4x4, const Matrix4x4 & viewportMatrix, uint32_t color);
	void DrawGrid(const Matrix4x4 & viewProjectionMatrix, const Matrix4x4 & viewportMatrix);
	void DrawPlane(const Plane & plane, const Matrix4x4 & viewProjectionMatrix, const Matrix4x4 & viewportMatrix, uint32_t color);
	void DrawRay(const Line & ray, float length, const Matrix4x4 & viewProjectionMatrix, const Matrix4x4 & viewportMatrix, uint32_t color);
	void DrawTriangle(const Triangle & triangle, const Matrix4x4 & viewProjectionMatrix, const Matrix4x4 & viewportMatrix, uint32_t color);

	Vector3 project = Project(Subtract(point, segment.origin), segment.diff);
	Vector3 closestPoint = ClosestPoint(point, segment);
	Sphere pointSphere{ point, 1.0f };
	Sphere closestPointSphere{ closestPoint, 1.0f };
	Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, 1280.0f, 720.0f, 0.0f, 1.0f);
	Matrix4x4 worldMatrix = MakeAffineMatrix({ 1.0f, 1.0f, 1.0f }, rotate, translate);
	Matrix4x4 cameraMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, { 0.0f,0.0f,0.0f });
	Matrix4x4 viewMatrix = Inverse(cameraMatrix);
	Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, 1.3f, 0.1f, 100.0f);
	Matrix4x4 viewProjectionMatrix = Multiply(worldMatrix, Multiply(viewMatrix, projectionMatrix));
	Vector3 start = Transform(Transform(segment.origin, viewProjectionMatrix), viewportMatrix);
	Vector3 end = Transform(Transform(Add(segment.origin, segment.diff), viewProjectionMatrix), viewportMatrix);
	plane.normal = Normalize(plane.normal);


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

		ImGui::Begin("Window");
		ImGui::DragFloat3("CameraTranslate", &cameraTranslate.x, 0.01f);
		ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.01f);
		ImGui::DragFloat3("SphereCenter", &pointSphere.center.x, 0.01f);
		ImGui::DragFloat("SphereRadius", &pointSphere.radius, 0.01f);
		ImGui::DragFloat3("planeNormal", &plane.normal.x, 0.01f);
		ImGui::DragFloat3("Segment", &line.origin.x, 0.01f);
		ImGui::DragFloat("planeDistance", &plane.distance, 0.01f);
		ImGui::End();

		worldMatrix = MakeAffineMatrix({ 1.0f, 1.0f, 1.0f }, rotate, translate);
		cameraMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, cameraTranslate);
		viewMatrix = Inverse(cameraMatrix);
		projectionMatrix = MakePerspectiveFovMatrix(0.45f, 1.3f, 0.1f, 100.0f);
		viewProjectionMatrix = Multiply(worldMatrix, Multiply(viewMatrix, projectionMatrix));

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		/*Novice::DrawLine(static_cast<int>(start.x), static_cast<int>(start.y),
			static_cast<int>(end.x), static_cast<int>(end.y), WHITE);*/
		DrawGrid(viewProjectionMatrix, viewportMatrix);
		if (IsCollision(triangle, segment)) {
			//DrawSphere(pointSphere, viewProjectionMatrix, viewportMatrix, RED);
			//DrawPlane(plane, viewProjectionMatrix, viewportMatrix, RED);
			DrawRay(line, 2.0f, viewProjectionMatrix, viewportMatrix, RED);
			DrawTriangle(triangle, viewProjectionMatrix, viewportMatrix, RED);
			//DrawSphere(closestPointSphere, viewProjectionMatrix, viewportMatrix, RED);

		} else {
			//DrawSphere(pointSphere, viewProjectionMatrix, viewportMatrix, WHITE);
			//DrawPlane(plane, viewProjectionMatrix, viewportMatrix, WHITE);
			DrawRay(line, 2.0f, viewProjectionMatrix, viewportMatrix, WHITE);
			DrawTriangle(triangle, viewProjectionMatrix, viewportMatrix, WHITE);
			//DrawSphere(closestPointSphere, viewProjectionMatrix, viewportMatrix, WHITE);
		}

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

Vector3 Add(const Vector3& v1, const Vector3& v2) {
	Vector3 result;

	result.x = v1.x + v2.x;
	result.y = v1.y + v2.y;
	result.z = v1.z + v2.z;

	return result;
}

Vector3 Subtract(const Vector3& v1, const Vector3& v2) {
	Vector3 result;

	result.x = v1.x - v2.x;
	result.y = v1.y - v2.y;
	result.z = v1.z - v2.z;

	return result;
}

Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result;

	for (int row = 0; row < 4; ++row) {
		for (int col = 0; col < 4; ++col) {
			result.m[row][col] = 0.0f;
			for (int k = 0; k < 4; ++k) {
				result.m[row][col] += m1.m[row][k] * m2.m[k][col];
			}
		}
	}

	return result;
}

Vector3 Multiply(float scalar, const Vector3& v1) {
	Vector3 result;

	result.x = v1.x * scalar;
	result.y = v1.y * scalar;
	result.z = v1.z * scalar;

	return result;
}

Matrix4x4 Inverse(const Matrix4x4& m1) {
	Matrix4x4 result = {};
	float det = 0.0f;

	// 行列式を計算
	det = m1.m[0][0] * (m1.m[1][1] * (m1.m[2][2] * m1.m[3][3] - m1.m[2][3] * m1.m[3][2]) -
		m1.m[1][2] * (m1.m[2][1] * m1.m[3][3] - m1.m[2][3] * m1.m[3][1]) +
		m1.m[1][3] * (m1.m[2][1] * m1.m[3][2] - m1.m[2][2] * m1.m[3][1])) -
		m1.m[0][1] * (m1.m[1][0] * (m1.m[2][2] * m1.m[3][3] - m1.m[2][3] * m1.m[3][2]) -
			m1.m[1][2] * (m1.m[2][0] * m1.m[3][3] - m1.m[2][3] * m1.m[3][0]) +
			m1.m[1][3] * (m1.m[2][0] * m1.m[3][2] - m1.m[2][2] * m1.m[3][0])) +
		m1.m[0][2] * (m1.m[1][0] * (m1.m[2][1] * m1.m[3][3] - m1.m[2][3] * m1.m[3][1]) -
			m1.m[1][1] * (m1.m[2][0] * m1.m[3][3] - m1.m[2][3] * m1.m[3][0]) +
			m1.m[1][3] * (m1.m[2][0] * m1.m[3][1] - m1.m[2][1] * m1.m[3][0])) -
		m1.m[0][3] * (m1.m[1][0] * (m1.m[2][1] * m1.m[3][2] - m1.m[2][2] * m1.m[3][1]) -
			m1.m[1][1] * (m1.m[2][0] * m1.m[3][2] - m1.m[2][2] * m1.m[3][0]) +
			m1.m[1][2] * (m1.m[2][0] * m1.m[3][1] - m1.m[2][1] * m1.m[3][0]));

	// 逆行列を計算
	float invDet = 1.0f / det;

	result.m[0][0] = invDet * (m1.m[1][1] * (m1.m[2][2] * m1.m[3][3] - m1.m[2][3] * m1.m[3][2]) -
		m1.m[1][2] * (m1.m[2][1] * m1.m[3][3] - m1.m[2][3] * m1.m[3][1]) +
		m1.m[1][3] * (m1.m[2][1] * m1.m[3][2] - m1.m[2][2] * m1.m[3][1]));
	result.m[0][1] = -invDet * (m1.m[0][1] * (m1.m[2][2] * m1.m[3][3] - m1.m[2][3] * m1.m[3][2]) -
		m1.m[0][2] * (m1.m[2][1] * m1.m[3][3] - m1.m[2][3] * m1.m[3][1]) +
		m1.m[0][3] * (m1.m[2][1] * m1.m[3][2] - m1.m[2][2] * m1.m[3][1]));
	result.m[0][2] = invDet * (m1.m[0][1] * (m1.m[1][2] * m1.m[3][3] - m1.m[1][3] * m1.m[3][2]) -
		m1.m[0][2] * (m1.m[1][1] * m1.m[3][3] - m1.m[1][3] * m1.m[3][1]) +
		m1.m[0][3] * (m1.m[1][1] * m1.m[3][2] - m1.m[1][2] * m1.m[3][1]));
	result.m[0][3] = -invDet * (m1.m[0][1] * (m1.m[1][2] * m1.m[2][3] - m1.m[1][3] * m1.m[2][2]) -
		m1.m[0][2] * (m1.m[1][1] * m1.m[2][3] - m1.m[1][3] * m1.m[2][1]) +
		m1.m[0][3] * (m1.m[1][1] * m1.m[2][2] - m1.m[1][2] * m1.m[2][1]));

	result.m[1][0] = -invDet * (m1.m[1][0] * (m1.m[2][2] * m1.m[3][3] - m1.m[2][3] * m1.m[3][2]) -
		m1.m[1][2] * (m1.m[2][0] * m1.m[3][3] - m1.m[2][3] * m1.m[3][0]) +
		m1.m[1][3] * (m1.m[2][0] * m1.m[3][2] - m1.m[2][2] * m1.m[3][0]));
	result.m[1][1] = invDet * (m1.m[0][0] * (m1.m[2][2] * m1.m[3][3] - m1.m[2][3] * m1.m[3][2]) -
		m1.m[0][2] * (m1.m[2][0] * m1.m[3][3] - m1.m[2][3] * m1.m[3][0]) +
		m1.m[0][3] * (m1.m[2][0] * m1.m[3][2] - m1.m[2][2] * m1.m[3][0]));
	result.m[1][2] = -invDet * (m1.m[0][0] * (m1.m[1][2] * m1.m[3][3] - m1.m[1][3] * m1.m[3][2]) -
		m1.m[0][2] * (m1.m[1][0] * m1.m[3][3] - m1.m[1][3] * m1.m[3][0]) +
		m1.m[0][3] * (m1.m[1][0] * m1.m[3][2] - m1.m[1][2] * m1.m[3][0]));
	result.m[1][3] = invDet * (m1.m[0][0] * (m1.m[1][2] * m1.m[2][3] - m1.m[1][3] * m1.m[2][2]) -
		m1.m[0][2] * (m1.m[1][0] * m1.m[2][3] - m1.m[1][3] * m1.m[2][0]) +
		m1.m[0][3] * (m1.m[1][0] * m1.m[2][2] - m1.m[1][2] * m1.m[2][0]));

	result.m[2][0] = invDet * (m1.m[1][0] * (m1.m[2][1] * m1.m[3][3] - m1.m[2][3] * m1.m[3][2]) -
		m1.m[1][1] * (m1.m[2][0] * m1.m[3][3] - m1.m[2][3] * m1.m[3][0]) +
		m1.m[1][3] * (m1.m[2][0] * m1.m[3][1] - m1.m[2][1] * m1.m[3][0]));
	result.m[2][1] = -invDet * (m1.m[0][0] * (m1.m[2][1] * m1.m[3][3] - m1.m[2][3] * m1.m[3][2]) -
		m1.m[0][1] * (m1.m[2][0] * m1.m[3][3] - m1.m[2][3] * m1.m[3][0]) +
		m1.m[0][3] * (m1.m[2][0] * m1.m[3][1] - m1.m[2][1] * m1.m[3][0]));
	result.m[2][2] = invDet * (m1.m[0][0] * (m1.m[1][1] * m1.m[3][3] - m1.m[1][3] * m1.m[3][1]) -
		m1.m[0][1] * (m1.m[1][0] * m1.m[3][3] - m1.m[1][3] * m1.m[3][0]) +
		m1.m[0][3] * (m1.m[1][0] * m1.m[3][1] - m1.m[1][1] * m1.m[3][0]));
	result.m[2][3] = -invDet * (m1.m[0][0] * (m1.m[1][1] * m1.m[2][3] - m1.m[1][3] * m1.m[2][1]) -
		m1.m[0][1] * (m1.m[1][0] * m1.m[2][3] - m1.m[1][3] * m1.m[2][0]) +
		m1.m[0][3] * (m1.m[1][0] * m1.m[2][1] - m1.m[1][1] * m1.m[2][0]));

	result.m[3][0] = -invDet * (m1.m[1][0] * (m1.m[2][1] * m1.m[3][2] - m1.m[2][2] * m1.m[3][1]) -
		m1.m[1][1] * (m1.m[2][0] * m1.m[3][2] - m1.m[2][2] * m1.m[3][0]) +
		m1.m[1][2] * (m1.m[2][0] * m1.m[3][1] - m1.m[2][1] * m1.m[3][0]));
	result.m[3][1] = invDet * (m1.m[0][0] * (m1.m[2][1] * m1.m[3][2] - m1.m[2][2] * m1.m[3][1]) -
		m1.m[0][1] * (m1.m[2][0] * m1.m[3][2] - m1.m[2][2] * m1.m[3][0]) +
		m1.m[0][2] * (m1.m[2][0] * m1.m[3][1] - m1.m[2][1] * m1.m[3][0]));
	result.m[3][2] = -invDet * (m1.m[0][0] * (m1.m[1][1] * m1.m[3][2] - m1.m[1][2] * m1.m[3][1]) -
		m1.m[0][1] * (m1.m[1][0] * m1.m[3][2] - m1.m[1][2] * m1.m[3][0]) +
		m1.m[0][2] * (m1.m[1][0] * m1.m[3][1] - m1.m[1][1] * m1.m[3][0]));
	result.m[3][3] = invDet * (m1.m[0][0] * (m1.m[1][1] * m1.m[2][2] - m1.m[1][2] * m1.m[2][1]) -
		m1.m[0][1] * (m1.m[1][0] * m1.m[2][2] - m1.m[1][2] * m1.m[2][0]) +
		m1.m[0][2] * (m1.m[1][0] * m1.m[2][1] - m1.m[1][1] * m1.m[2][0]));

	return result;
}

float Dot(const Vector3& v1, const Vector3& v2) {
	float result;

	result = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;

	return result;
}

Vector3 Normalize(Vector3 v1) {
	Vector3 result;

	result.x = v1.x / std::sqrtf(Dot(v1, v1));
	result.y = v1.y / std::sqrtf(Dot(v1, v1));
	result.z = v1.z / std::sqrtf(Dot(v1, v1));

	return result;
}

Vector3 Cross(const Vector3& v1, const Vector3& v2) {
	Vector3 result;
	result.x = v1.y * v2.z - v1.z * v2.y;
	result.y = v1.z * v2.x - v1.x * v2.z;
	result.z = v1.x * v2.y - v1.y * v2.x;
	return result;
}

Vector3 Project(const Vector3& v1, const Vector3& v2) {
	float v2LengthSq = v2.x * v2.x + v2.y * v2.y + v2.z * v2.z;
	assert(v2LengthSq != 0.0f);
	float dot = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	float t = dot / v2LengthSq;
	return { v2.x * t, v2.y * t, v2.z * t };
}

Vector3 ClosestPoint(const Vector3& point, const Segment& segment) {
	// 線分の始点と方向
	Vector3 segOrigin = segment.origin;
	Vector3 segDiff = segment.diff;

	// 始点から点へのベクトル
	Vector3 v = Subtract(point, segOrigin);

	// 線分方向ベクトルの長さの2乗
	float dLenSq = segDiff.x * segDiff.x + segDiff.y * segDiff.y + segDiff.z * segDiff.z;
	if (dLenSq == 0.0f) {
		// 長さ0の線分の場合は始点を返す
		return segOrigin;
	}

	// 内積
	float t = (v.x * segDiff.x + v.y * segDiff.y + v.z * segDiff.z) / dLenSq;

	// 0～1にクランプ
	t = std::clamp(t, 0.0f, 1.0f);

	// 最近接点
	Vector3 result;
	result.x = segOrigin.x + segDiff.x * t;
	result.y = segOrigin.y + segDiff.y * t;
	result.z = segOrigin.z + segDiff.z * t;
	return result;
}

Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth) {
	Matrix4x4 result;

	result.m[0][0] = width / 2.0f;
	result.m[0][1] = 0.0f;
	result.m[0][2] = 0.0f;
	result.m[0][3] = 0.0f;

	result.m[1][0] = 0.0f;
	result.m[1][1] = -(height / 2.0f);
	result.m[1][2] = 0.0f;
	result.m[1][3] = 0.0f;

	result.m[2][0] = 0.0f;
	result.m[2][1] = 0.0f;
	result.m[2][2] = maxDepth - minDepth;
	result.m[2][3] = 0.0f;

	result.m[3][0] = left + width / 2.0f;
	result.m[3][1] = top + height / 2.0f;
	result.m[3][2] = minDepth;
	result.m[3][3] = 1.0f;

	return result;
}

Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {
	Matrix4x4 result = {};

	float yScale = 1.0f / std::tan(fovY / 2.0f);
	float xScale = yScale / aspectRatio;

	result.m[0][0] = xScale;
	result.m[1][1] = yScale;
	result.m[2][2] = farClip / (farClip - nearClip);
	result.m[2][3] = 1.0f;
	result.m[3][2] = -nearClip * farClip / (farClip - nearClip);

	return result;
}

Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {

	Matrix4x4 MakeRotateXMatrix;

	MakeRotateXMatrix.m[0][0] = 1.0f;
	MakeRotateXMatrix.m[0][1] = 0.0f;
	MakeRotateXMatrix.m[0][2] = 0.0f;
	MakeRotateXMatrix.m[0][3] = 0.0f;

	MakeRotateXMatrix.m[1][0] = 0.0f;
	MakeRotateXMatrix.m[1][1] = std::cos(rotate.x);
	MakeRotateXMatrix.m[1][2] = std::sin(rotate.x);
	MakeRotateXMatrix.m[1][3] = 0.0f;

	MakeRotateXMatrix.m[2][0] = 0.0f;
	MakeRotateXMatrix.m[2][1] = -std::sin(rotate.x);
	MakeRotateXMatrix.m[2][2] = std::cos(rotate.x);
	MakeRotateXMatrix.m[2][3] = 0.0f;

	MakeRotateXMatrix.m[3][0] = 0.0f;
	MakeRotateXMatrix.m[3][1] = 0.0f;
	MakeRotateXMatrix.m[3][2] = 0.0f;
	MakeRotateXMatrix.m[3][3] = 1.0f;

	Matrix4x4 MakeRotateYMatrix;

	MakeRotateYMatrix.m[0][0] = std::cos(rotate.y);
	MakeRotateYMatrix.m[0][1] = 0.0f;
	MakeRotateYMatrix.m[0][2] = -std::sin(rotate.y);
	MakeRotateYMatrix.m[0][3] = 0.0f;

	MakeRotateYMatrix.m[1][0] = 0.0f;
	MakeRotateYMatrix.m[1][1] = 1.0f;
	MakeRotateYMatrix.m[1][2] = 0.0f;
	MakeRotateYMatrix.m[1][3] = 0.0f;

	MakeRotateYMatrix.m[2][0] = std::sin(rotate.y);
	MakeRotateYMatrix.m[2][1] = 0.0f;
	MakeRotateYMatrix.m[2][2] = std::cos(rotate.y);
	MakeRotateYMatrix.m[2][3] = 0.0f;

	MakeRotateYMatrix.m[3][0] = 0.0f;
	MakeRotateYMatrix.m[3][1] = 0.0f;
	MakeRotateYMatrix.m[3][2] = 0.0f;
	MakeRotateYMatrix.m[3][3] = 1.0f;

	Matrix4x4 MakeRotateZMatrix;

	MakeRotateZMatrix.m[0][0] = std::cos(rotate.z);
	MakeRotateZMatrix.m[0][1] = std::sin(rotate.z);
	MakeRotateZMatrix.m[0][2] = 0.0f;
	MakeRotateZMatrix.m[0][3] = 0.0f;

	MakeRotateZMatrix.m[1][0] = -std::sin(rotate.z);
	MakeRotateZMatrix.m[1][1] = std::cos(rotate.z);
	MakeRotateZMatrix.m[1][2] = 0.0f;
	MakeRotateZMatrix.m[1][3] = 0.0f;

	MakeRotateZMatrix.m[2][0] = 0.0f;
	MakeRotateZMatrix.m[2][1] = 0.0f;
	MakeRotateZMatrix.m[2][2] = 1.0f;
	MakeRotateZMatrix.m[2][3] = 0.0f;

	MakeRotateZMatrix.m[3][0] = 0.0f;
	MakeRotateZMatrix.m[3][1] = 0.0f;
	MakeRotateZMatrix.m[3][2] = 0.0f;
	MakeRotateZMatrix.m[3][3] = 1.0f;

	Matrix4x4 MakeRotateXYZMatrix = Multiply(MakeRotateXMatrix, Multiply(MakeRotateYMatrix, MakeRotateZMatrix));

	Matrix4x4 result;

	result.m[0][0] = MakeRotateXYZMatrix.m[0][0] * scale.x;
	result.m[0][1] = MakeRotateXYZMatrix.m[0][1] * scale.x;
	result.m[0][2] = MakeRotateXYZMatrix.m[0][2] * scale.x;
	result.m[0][3] = 0.0f;

	result.m[1][0] = MakeRotateXYZMatrix.m[1][0] * scale.y;
	result.m[1][1] = MakeRotateXYZMatrix.m[1][1] * scale.y;
	result.m[1][2] = MakeRotateXYZMatrix.m[1][2] * scale.y;
	result.m[1][3] = 0.0f;

	result.m[2][0] = MakeRotateXYZMatrix.m[2][0] * scale.z;
	result.m[2][1] = MakeRotateXYZMatrix.m[2][1] * scale.z;
	result.m[2][2] = MakeRotateXYZMatrix.m[2][2] * scale.z;
	result.m[2][3] = 0.0f;

	result.m[3][0] = translate.x;
	result.m[3][1] = translate.y;
	result.m[3][2] = translate.z;
	result.m[3][3] = 1.0f;

	return result;
}

Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {
	Vector3 result;

	result.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + 1.0f * matrix.m[3][0];
	result.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + 1.0f * matrix.m[3][1];
	result.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + 1.0f * matrix.m[3][2];
	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + 1.0f * matrix.m[3][3];


	assert(w != 0.0f);
	result.x /= w;
	result.y /= w;
	result.z /= w;

	return result;
}

bool IsCollision(const Sphere& s1, const Sphere& s2) {
	float dx = s1.center.x - s2.center.x;
	float dy = s1.center.y - s2.center.y;
	float dz = s1.center.z - s2.center.z;
	float distanceSq = dx * dx + dy * dy + dz * dz;
	float radiusSum = s1.radius + s2.radius;
	return distanceSq <= radiusSum * radiusSum;
}

bool IsCollision(const Sphere& s, const Plane& p) {
	float distance = Dot(s.center, p.normal) - p.distance;
	return distance <= s.radius;
}

bool IsCollision(const Segment& segment, const Plane& plane) {
	// まず垂直判定を行うために法線と線の内積を求める
	float dot = Dot(plane.normal, segment.diff);
	if (dot == 0.0f) {
		// 線分が平面と平行
		return false;
	}
	// tを求める
	float t = (plane.distance - Dot(segment.origin, plane.normal)) / dot;

	// tの値が0～1の範囲なら線分と平面は交差している
	return (0.0f <= t && t <= 1.0f);
}

bool IsCollision(const Triangle& triangle, const Segment& segment) {
	// 三角形の3頂点
	const Vector3& v0 = triangle.vertices[0];
	const Vector3& v1 = triangle.vertices[1];
	const Vector3& v2 = triangle.vertices[2];

	// 三角形の法線を計算
	Vector3 edge1 = Subtract(v1, v0);
	Vector3 edge2 = Subtract(v2, v0);
	Vector3 normal = Normalize(Cross(edge1, edge2));

	// 平面の方程式: normal・(P - v0) = 0
	// 線分の始点と終点
	Vector3 p0 = segment.origin;
	Vector3 p1 = Add(segment.origin, segment.diff);

	// 線分の方向
	Vector3 dir = Subtract(p1, p0);

	float denom = Dot(normal, dir);
	if (std::fabs(denom) < 1e-6f) {
		// 線分が平面と平行
		return false;
	}

	float t = Dot(normal, Subtract(v0, p0)) / denom;
	if (t < 0.0f || t > 1.0f) {
		// 交点が線分の範囲外
		return false;
	}

	// 交点を求める
	Vector3 intersect = Add(p0, Multiply(t, dir));

	// バリツェ法で三角形内か判定
	Vector3 c0 = Cross(Subtract(v1, v0), Subtract(intersect, v0));
	Vector3 c1 = Cross(Subtract(v2, v1), Subtract(intersect, v1));
	Vector3 c2 = Cross(Subtract(v0, v2), Subtract(intersect, v2));

	if (Dot(normal, c0) >= 0 && Dot(normal, c1) >= 0 && Dot(normal, c2) >= 0) {
		return true;
	}
	return false;
}

Vector3 Perpendicular(const Vector3& vector) {
	if (vector.x != 0.0f || vector.y != 0.0f) {
		return { -vector.y, vector.x, 0.0f };
	}
	return { 0.0f, -vector.z, vector.y };
}

void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	const uint32_t kSubdivision = 18;   // 分割数（緯度・経度）
	const float kLonEvery = 2.0f * 3.14159265f / kSubdivision; // 経度分割単位角（0〜2π）
	const float kLatEvery = 3.14159265f / kSubdivision;        // 緯度分割単位角（-π/2〜π/2）

	for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex) {
		float lat = -3.14159265f / 2.0f + kLatEvery * latIndex;           // 現在の緯度
		float nextLat = lat + kLatEvery;

		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex) {
			float lon = kLonEvery * lonIndex;                             // 現在の経度
			float nextLon = lon + kLonEvery;

			// ワールド座標の3点（a: 緯度×経度, b: 緯度×次経度, c: 次緯度×経度）
			Vector3 a = {
				sphere.center.x + sphere.radius * cosf(lat) * cosf(lon),
				sphere.center.y + sphere.radius * sinf(lat),
				sphere.center.z + sphere.radius * cosf(lat) * sinf(lon)
			};
			Vector3 b = {
				sphere.center.x + sphere.radius * cosf(lat) * cosf(nextLon),
				sphere.center.y + sphere.radius * sinf(lat),
				sphere.center.z + sphere.radius * cosf(lat) * sinf(nextLon)
			};
			Vector3 c = {
				sphere.center.x + sphere.radius * cosf(nextLat) * cosf(lon),
				sphere.center.y + sphere.radius * sinf(nextLat),
				sphere.center.z + sphere.radius * cosf(nextLat) * sinf(lon)
			};

			// 画面座標に変換
			Vector3 screenA = Transform(Transform(a, viewProjectionMatrix), viewportMatrix);
			Vector3 screenB = Transform(Transform(b, viewProjectionMatrix), viewportMatrix);
			Vector3 screenC = Transform(Transform(c, viewProjectionMatrix), viewportMatrix);

			// 線を描画（ab 緯度方向, ac 経度方向）
			Novice::DrawLine((int)screenA.x, (int)screenA.y, (int)screenB.x, (int)screenB.y, color);
			Novice::DrawLine((int)screenA.x, (int)screenA.y, (int)screenC.x, (int)screenC.y, color);
		}
	}
}

void DrawPlane(const Plane& plane, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	Vector3 center = Multiply(plane.distance, plane.normal);
	Vector3 perpendiculars[4];
	perpendiculars[0] = Normalize(Perpendicular(plane.normal));
	perpendiculars[1] = { -perpendiculars[0].x, -perpendiculars[0].y, -perpendiculars[0].z };
	perpendiculars[2] = Cross(plane.normal, perpendiculars[0]);
	perpendiculars[3] = { -perpendiculars[2].x, -perpendiculars[2].y, -perpendiculars[2].z };

	Vector3 points[4];
	for (int32_t index = 0; index < 4; ++index) {
		Vector3 extend = Multiply(2.0f, perpendiculars[index]);
		Vector3 point = Add(center, extend);
		points[index] = Transform(Transform(point, viewProjectionMatrix), viewportMatrix);
	}
	//pointsを結んでDrawLineで矩形を描画
	for (int32_t index = 0; index < 4; ++index) {
		Vector3 start = points[index];
		Vector3 end = points[(index + 1) % 4];
		Novice::DrawLine(static_cast<int>(start.x), static_cast<int>(start.y),
			static_cast<int>(end.x), static_cast<int>(end.y), color);
	}

}

void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {
	const float kGridHalfWidth = 2.0f; // グリッドの半幅（-2.0f〜+2.0f）
	const uint32_t kSubdivision = 10;  // 分割数（11本の線）
	const float kGridEvery = (kGridHalfWidth * 2.0f) / static_cast<float>(kSubdivision);

	// X方向の線（Zが固定でXが変わる）
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {
		float x = -kGridHalfWidth + kGridEvery * xIndex;

		Vector3 worldStart = { x, 0.0f, -kGridHalfWidth };
		Vector3 worldEnd = { x, 0.0f, +kGridHalfWidth };

		Vector3 ndcStart = Transform(worldStart, viewProjectionMatrix);
		Vector3 ndcEnd = Transform(worldEnd, viewProjectionMatrix);

		Vector3 screenStart = Transform(ndcStart, viewportMatrix);
		Vector3 screenEnd = Transform(ndcEnd, viewportMatrix);

		Novice::DrawLine(static_cast<int>(screenStart.x), static_cast<int>(screenStart.y),
			static_cast<int>(screenEnd.x), static_cast<int>(screenEnd.y),
			BLACK);
	}

	// Z方向の線（Xが固定でZが変わる）
	for (uint32_t zIndex = 0; zIndex <= kSubdivision; ++zIndex) {
		float z = -kGridHalfWidth + kGridEvery * zIndex;

		Vector3 worldStart = { -kGridHalfWidth, 0.0f, z };
		Vector3 worldEnd = { +kGridHalfWidth, 0.0f, z };

		Vector3 ndcStart = Transform(worldStart, viewProjectionMatrix);
		Vector3 ndcEnd = Transform(worldEnd, viewProjectionMatrix);

		Vector3 screenStart = Transform(ndcStart, viewportMatrix);
		Vector3 screenEnd = Transform(ndcEnd, viewportMatrix);

		Novice::DrawLine(static_cast<int>(screenStart.x), static_cast<int>(screenStart.y),
			static_cast<int>(screenEnd.x), static_cast<int>(screenEnd.y),
			BLACK);
	}
}

void DrawRay(const Line& ray, float length, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	// 始点
	Vector3 start = Transform(Transform(ray.origin, viewProjectionMatrix), viewportMatrix);
	// 終点（始点 + 方向ベクトル * 長さ）
	Vector3 dir = Normalize(ray.diff);
	Vector3 endPoint3D = Add(ray.origin, Multiply(length, dir));
	Vector3 end = Transform(Transform(endPoint3D, viewProjectionMatrix), viewportMatrix);

	Novice::DrawLine(
		static_cast<int>(start.x), static_cast<int>(start.y),
		static_cast<int>(end.x), static_cast<int>(end.y),
		color
	);
}

void DrawTriangle(const Triangle& triangle, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	// 3頂点をワールド→ビュー→スクリーン座標へ変換
	Vector3 screenVertices[3];
	for (int i = 0; i < 3; ++i) {
		screenVertices[i] = Transform(Transform(triangle.vertices[i], viewProjectionMatrix), viewportMatrix);
	}

	// 3辺を線で描画
	for (int i = 0; i < 3; ++i) {
		const Vector3& v0 = screenVertices[i];
		const Vector3& v1 = screenVertices[(i + 1) % 3];
		Novice::DrawLine(
			static_cast<int>(v0.x), static_cast<int>(v0.y),
			static_cast<int>(v1.x), static_cast<int>(v1.y),
			color
		);
	}

}