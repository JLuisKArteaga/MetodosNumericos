#include "raylib.h"
#include "raymath.h"
#include <cmath>
#include <vector>
#include <deque>

// ============================================
// CONFIGURACION
// ============================================
const int MAX_PARTICULAS = 500;
const int MAX_TRAIL = 50;

float sigma = 10.0f;
float rho = 28.0f;
float beta = 8.0f / 3.0f;
float dt = 0.01f;
float factorLorenz = 0.04f;
float fuerzaSpin = 8.0f;

// ============================================
// CÁMARA MANUAL
// ============================================
struct CamaraControl {
    float rotX = 30.0f;
    float rotY = 45.0f;
    float distancia = 70.0f;
    Vector3 objetivo = { 0, 0, 25 };
    bool arrastrando = false;
    Vector2 posAnterior;

    void actualizar() {
        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            arrastrando = true;
            posAnterior = GetMousePosition();
        }
        if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) {
            arrastrando = false;
        }

        if (arrastrando) {
            Vector2 mouse = GetMousePosition();
            rotY += (mouse.x - posAnterior.x) * 0.5f;
            rotX -= (mouse.y - posAnterior.y) * 0.5f;

            if (rotX > 85.0f) rotX = 85.0f;
            if (rotX < -85.0f) rotX = -85.0f;
            posAnterior = mouse;
        }

        float rueda = GetMouseWheelMove();
        distancia -= rueda * 3.0f;
        distancia = fmaxf(20.0f, fminf(150.0f, distancia));
    }

    Camera3D obtenerCamara() {
        float rx = rotX * DEG2RAD;
        float ry = rotY * DEG2RAD;

        Camera3D cam = { 0 };
        cam.position = {
            objetivo.x + distancia * cosf(rx) * sinf(ry),
            objetivo.y + distancia * sinf(rx),
            objetivo.z + distancia * cosf(rx) * cosf(ry)
        };
        cam.target = objetivo;
        cam.up = { 0, 1, 0 };
        cam.fovy = 60.0f;
        cam.projection = CAMERA_PERSPECTIVE;
        return cam;
    }
};

// ============================================
// PARTÍCULA CON RK4 + TORQUE DIFERENCIAL
// ============================================
struct Particula {
    Vector3 pos;
    Vector3 vel;
    std::deque<Vector3> rastro;
    Color color;
    float vida;
    bool activa;
    int tipo; // 0 = flujo a agujero, 1 = órbita planeta

    Particula() : vida(0), activa(false), tipo(0) {}

    void reset(Vector3 p, Vector3 v, Color c, int t) {
        pos = p;
        vel = v;
        color = c;
        vida = 0;
        activa = true;
        tipo = t;
        rastro.clear();
        rastro.push_back(p);
    }

    Vector3 lorenz(Vector3 p) {
        return {
            sigma * (p.y - p.x),
            p.x * (rho - p.z) - p.y,
            p.x * p.y - beta * p.z
        };
    }

    void calcularFuerzas(Vector3 p, Vector3& fLorenz, Vector3& fGrav,
        Vector3& fSpin, Vector3 planeta, Vector3 agujero) {

        fLorenz = lorenz(p);

        Vector3 haciaAgujero = Vector3Subtract(agujero, p);
        float distAgujero = Vector3Length(haciaAgujero);
        Vector3 dirAgujero = Vector3Normalize(haciaAgujero);

        Vector3 haciaPlaneta = Vector3Subtract(planeta, p);
        float distPlaneta = Vector3Length(haciaPlaneta);

        fGrav = { 0, 0, 0 };

        // Atracción al agujero (más débil para tipo 1)
        float masaAgujero = (tipo == 0) ? 500.0f : 150.0f;
        if (distAgujero > 1.0f) {
            fGrav = Vector3Add(fGrav,
                Vector3Scale(dirAgujero, masaAgujero / (distAgujero * distAgujero + 10.0f)));
        }

        // Atracción al planeta (más fuerte para tipo 1)
        float masaPlaneta = (tipo == 0) ? 50.0f : 200.0f;
        if (distPlaneta > 0.5f) {
            Vector3 dirPlaneta = Vector3Normalize(haciaPlaneta);
            fGrav = Vector3Add(fGrav,
                Vector3Scale(dirPlaneta, masaPlaneta / (distPlaneta * distPlaneta)));
        }

        // TORQUE DIFERENCIAL
        fSpin = { 0, 0, 0 };

        // Rotación ANTI-HORARIA cerca del planeta (más fuerte para tipo 1)
        if (distPlaneta < 25.0f) {
            Vector3 spinA = { -haciaPlaneta.y, haciaPlaneta.x, 0 };
            float factorA = 1.0f - (distPlaneta / 25.0f);
            float intensidad = (tipo == 0) ? fuerzaSpin : fuerzaSpin * 2.5f;
            fSpin = Vector3Add(fSpin,
                Vector3Scale(Vector3Normalize(spinA), intensidad * factorA));
        }

        // Rotación HORARIA cerca del agujero
        if (distAgujero < 20.0f) {
            Vector3 spinB = { haciaAgujero.y, -haciaAgujero.x, 0 };
            float factorB = 1.0f - (distAgujero / 20.0f);
            fSpin = Vector3Add(fSpin,
                Vector3Scale(Vector3Normalize(spinB), fuerzaSpin * factorB));
        }
    }

    void actualizar(Vector3 planeta, Vector3 agujero) {
        if (!activa) return;

        auto evaluar = [&](Vector3 p, Vector3 v, Vector3& dvdt) {
            Vector3 fLorenz, fGrav, fSpin;
            calcularFuerzas(p, fLorenz, fGrav, fSpin, planeta, agujero);

            Vector3 acc = Vector3Add(
                Vector3Scale(fLorenz, factorLorenz),
                Vector3Add(fGrav, fSpin)
            );
            dvdt = acc;
            return v;
            };

        Vector3 k1_p, k2_p, k3_p, k4_p;
        Vector3 k1_v, k2_v, k3_v, k4_v;
        Vector3 p_temp, v_temp;

        k1_p = evaluar(pos, vel, k1_v);

        p_temp = Vector3Add(pos, Vector3Scale(k1_p, dt * 0.5f));
        v_temp = Vector3Add(vel, Vector3Scale(k1_v, dt * 0.5f));
        k2_p = evaluar(p_temp, v_temp, k2_v);

        p_temp = Vector3Add(pos, Vector3Scale(k2_p, dt * 0.5f));
        v_temp = Vector3Add(vel, Vector3Scale(k2_v, dt * 0.5f));
        k3_p = evaluar(p_temp, v_temp, k3_v);

        p_temp = Vector3Add(pos, Vector3Scale(k3_p, dt));
        v_temp = Vector3Add(vel, Vector3Scale(k3_v, dt));
        k4_p = evaluar(p_temp, v_temp, k4_v);

        Vector3 suma_p = Vector3Add(
            Vector3Add(k1_p, Vector3Scale(k2_p, 2.0f)),
            Vector3Add(Vector3Scale(k3_p, 2.0f), k4_p)
        );
        Vector3 suma_v = Vector3Add(
            Vector3Add(k1_v, Vector3Scale(k2_v, 2.0f)),
            Vector3Add(Vector3Scale(k3_v, 2.0f), k4_v)
        );

        pos = Vector3Add(pos, Vector3Scale(suma_p, dt / 6.0f));
        vel = Vector3Add(vel, Vector3Scale(suma_v, dt / 6.0f));

        vel = Vector3Scale(vel, 0.99f);
        vida += dt;

        rastro.push_front(pos);
        if (rastro.size() > MAX_TRAIL) rastro.pop_back();
    }

    bool debeReset(Vector3 planeta, Vector3 agujero) {
        if (!activa) return true;
        float da = Vector3Distance(pos, agujero);
        float dp = Vector3Distance(pos, planeta);
        float tiempoMax = (tipo == 0) ? 15.0f : 25.0f; // Órbitas duran más
        return (da < 0.8f) || (dp > 100.0f) || (vida > tiempoMax) || (da > 150.0f);
    }

    void dibujar(Vector3 agujero) {
        if (!activa || rastro.size() < 2) return;

        float d = Vector3Distance(pos, agujero);
        float t = fminf(d / 30.0f, 1.0f);

        // Color diferente según tipo
        Color c;
        if (tipo == 0) {
            // Naranja a azul (flujo normal)
            c = Color{
                (unsigned char)(255 * t + 255 * (1.0f - t)),
                (unsigned char)(150 * t + 100 * (1.0f - t)),
                (unsigned char)(255 * (1.0f - t)),
                255
            };
        }
        else {
            // Verde (órbita planeta)
            c = Color{
                (unsigned char)(100 * t),
                255,
                (unsigned char)(100 * (1.0f - t)),
                255
            };
        }

        for (size_t i = 0; i + 1 < rastro.size(); i++) {
            float alpha = 1.0f - ((float)i / rastro.size());
            Color cc = c;
            cc.a = (unsigned char)(255 * alpha);
            DrawLine3D(rastro[i], rastro[i + 1], cc);
        }

        DrawSphere(pos, (tipo == 0) ? 0.15f : 0.2f, WHITE);
    }
};

// ============================================
// SISTEMA CON DOS TIPOS DE PARTÍCULAS
// ============================================
struct Sistema {
    std::vector<Particula> particulas;
    Vector3 planeta;
    Vector3 agujero;

    void init() {
        particulas.resize(MAX_PARTICULAS);
        calcularFocos();
        for (int i = 0; i < MAX_PARTICULAS; i++) spawn(i);
    }

    void calcularFocos() {
        float b = sqrtf(fmaxf(0.0f, beta * (rho - 1.0f)));
        planeta = { -b, -b, rho - 1.0f };
        agujero = { b, b, rho - 1.0f };
    }

    void spawn(int idx) {
        // 20% de partículas orbitan el planeta, 80% fluyen al agujero
        int tipo = (GetRandomValue(0, 100) < 20) ? 1 : 0;

        float ang = GetRandomValue(0, 360) * DEG2RAD;
        Vector3 dir = Vector3Normalize(Vector3Subtract(agujero, planeta));
        Vector3 perp = { -dir.z, 0, dir.x };
        if (Vector3Length(perp) < 0.01f) perp = { 0, 1, 0 };
        perp = Vector3Normalize(perp);

        Vector3 perpRot = {
            perp.x * cosf(ang) - perp.z * sinf(ang),
            sinf(ang),
            perp.x * sinf(ang) + perp.z * cosf(ang)
        };

        Vector3 pos;
        Vector3 vel;
        Color c;

        if (tipo == 0) {
            // TIPO 0: Flujo hacia agujero (como antes)
            float radio = 1.2f + GetRandomValue(0, 300) / 100.0f;
            float altura = GetRandomValue(-150, 150) / 100.0f;

            pos = Vector3Add(planeta,
                Vector3Add(
                    Vector3Scale(dir, 1.0f + GetRandomValue(0, 100) / 100.0f),
                    Vector3Add(
                        Vector3Scale(perpRot, radio),
                        { 0, 0, altura }
                    )
                )
            );

            float velMag = 6.0f + GetRandomValue(0, 400) / 100.0f;
            Vector3 velDir = Vector3Normalize(Vector3Add(dir, Vector3Scale(perpRot, 0.3f)));
            vel = Vector3Scale(velDir, velMag);

            unsigned char intensidad = (unsigned char)(150 + velMag * 10);
            c = { 255, intensidad, (unsigned char)(50 + velMag * 20), 255 };
        }
        else {
            // TIPO 1: Órbita alrededor del planeta (anti-horario)
            float radio = 2.0f + GetRandomValue(0, 400) / 100.0f; // 2.0 a 6.0
            float altura = GetRandomValue(-100, 100) / 100.0f;

            // Posición en órbita circular
            pos = Vector3Add(planeta,
                Vector3Add(
                    Vector3Scale(perpRot, radio),
                    { 0, 0, altura }
                )
            );

            // Velocidad tangencial PURA (anti-horaria)
            // v = sqrt(GM/r) para órbita circular
            float velOrbita = sqrtf(200.0f / radio) * (0.8f + GetRandomValue(0, 40) / 100.0f);

            // Dirección tangencial anti-horaria
            Vector3 tangencial = { -perpRot.y, perpRot.x, 0 };
            tangencial = Vector3Normalize(tangencial);

            // Añadir pequeña componente hacia agujero para que eventualmente caigan
            vel = Vector3Add(
                Vector3Scale(tangencial, velOrbita),
                Vector3Scale(dir, 1.5f) // Lento acercamiento al agujero
            );

            c = { 100, 255, 100, 255 }; // Verde claro
        }

        particulas[idx].reset(pos, vel, c, tipo);
    }

    void update() {
        calcularFocos();
        for (int i = 0; i < MAX_PARTICULAS; i++) {
            particulas[i].actualizar(planeta, agujero);
            if (particulas[i].debeReset(planeta, agujero)) spawn(i);
        }
    }

    void dibujarCuerpos() {
        // Planeta
        DrawSphere(planeta, 1.5f, Color{ 255, 200, 50, 255 });
        DrawSphereWires(planeta, 1.6f, 16, 16, ORANGE);

        for (int i = 0; i < 4; i++) {
            float ang = i * 90 * DEG2RAD + GetTime() * 2.0f;
            Vector3 r = { cosf(ang) * 1.8f, sinf(ang) * 1.8f, 0 };
            DrawLine3D(planeta, Vector3Add(planeta, r), Color{ 255, 150, 0, 150 });
        }

        // Agujero negro
        for (int i = 0; i < 5; i++) {
            float r = 1.0f + i * 0.8f;
            float alpha = 1.0f - i * 0.15f;
            Color c = Color{
                (unsigned char)(255 * alpha),
                (unsigned char)(100 * alpha),
                (unsigned char)(255 * alpha),
                (unsigned char)(150 * alpha)
            };
            DrawCylinderEx(
                Vector3Add(agujero, { 0, 0, -0.05f * (5 - i) }),
                Vector3Add(agujero, { 0, 0, 0.05f * (5 - i) }),
                r, r, 32, c
            );
        }

        for (int i = 0; i < 4; i++) {
            float ang = -i * 90 * DEG2RAD - GetTime() * 3.0f;
            Vector3 r = { cosf(ang) * 2.2f, sinf(ang) * 2.2f, 0 };
            DrawLine3D(agujero, Vector3Add(agujero, r), Color{ 200, 0, 255, 150 });
        }

        DrawSphere(agujero, 0.8f, BLACK);
        DrawSphereWires(agujero, 0.9f, 16, 16, PURPLE);

        DrawCylinderEx(
            Vector3Add(agujero, { 0,0,1 }),
            Vector3Add(agujero, { 0,0,12 }),
            0.3f, 0, 8, Color{ 150, 0, 255, 100 }
        );
        DrawCylinderEx(
            Vector3Add(agujero, { 0,0,-1 }),
            Vector3Add(agujero, { 0,0,-12 }),
            0.3f, 0, 8, Color{ 150, 0, 255, 100 }
        );
    }

    void dibujarParticulas() {
        for (auto& p : particulas) p.dibujar(agujero);
    }
};

// ============================================
// SLIDER
// ============================================
struct Slider {
    Rectangle r;
    float minV, maxV, * val;
    const char* label;
    Color col;

    Slider(float x, float y, float w, float h, float mn, float mx, float* v, const char* l, Color c)
        : r{ x,y,w,h }, minV(mn), maxV(mx), val(v), label(l), col(c) {
    }

    void update() {
        if (IsMouseButtonDown(MOUSE_LEFT_BUTTON)) {
            Vector2 m = GetMousePosition();
            if (CheckCollisionPointRec(m, r)) {
                float pct = (m.x - r.x) / r.width;
                *val = minV + pct * (maxV - minV);
            }
        }
        *val = fmaxf(minV, fminf(maxV, *val));
    }

    void draw() {
        DrawRectangleRec(r, DARKGRAY);
        float pct = (*val - minV) / (maxV - minV);
        DrawRectangle(r.x, r.y, r.width * pct, r.height, col);
        DrawRectangleLinesEx(r, 2, WHITE);
        DrawText(label, (int)r.x, (int)r.y - 18, 16, RAYWHITE);
        DrawText(TextFormat("%.2f", *val), (int)(r.x + r.width + 5), (int)r.y, 14, YELLOW);
    }
};

// ============================================
// MAIN
// ============================================
int main() {
    InitWindow(1280, 720, "Lorenz - Doble Comportamiento");
    SetTargetFPS(60);

    CamaraControl cam;
    Sistema sys;
    sys.init();

    Slider s1(20, 60, 180, 18, 1, 30, &sigma, "Sigma", ORANGE);
    Slider s2(20, 110, 180, 18, 1, 50, &rho, "Rho", GREEN);
    Slider s3(20, 160, 180, 18, 0.1f, 10, &beta, "Beta", BLUE);
    Slider s4(20, 210, 180, 18, 0.001f, 0.05f, &dt, "dt", RED);
    Slider s5(20, 260, 180, 18, 0.0f, 0.5f, &factorLorenz, "Lorenz", PURPLE);
    Slider s6(20, 310, 180, 18, 0.0f, 20.0f, &fuerzaSpin, "Spin", YELLOW);
    Slider* sliders[] = { &s1, &s2, &s3, &s4, &s5, &s6 };

    bool pause = false;

    while (!WindowShouldClose()) {
        cam.actualizar();
        for (auto s : sliders) s->update();
        if (IsKeyPressed(KEY_SPACE)) pause = !pause;
        if (IsKeyPressed(KEY_R)) sys.init();
        if (!pause) sys.update();

        BeginDrawing();
        ClearBackground(Color{ 2, 2, 10, 255 });

        BeginMode3D(cam.obtenerCamara());
        DrawGrid(60, 2.0f);
        sys.dibujarCuerpos();
        sys.dibujarParticulas();
        EndMode3D();

        DrawRectangle(10, 10, 280, 360, Color{ 0,0,0,180 });
        DrawText("LORENZ - DOBLE COMPORTAMIENTO", 20, 15, 20, RAYWHITE);
        DrawText("Naranja: Flujo a agujero", 20, 38, 12, ORANGE);
        DrawText("Verde: Orbita planeta (anti-horario)", 20, 52, 12, GREEN);
        DrawText("Arrastra: Rotar | Rueda: Zoom", 20, 70, 11, GRAY);

        for (auto s : sliders) s->draw();

        DrawText(TextFormat("FPS: %d %s", GetFPS(), pause ? "[PAUSA]" : ""), 20, 340, 14, pause ? RED : GREEN);

        DrawCircle(300, 30, 8, Color{ 255, 200, 50, 255 });
        DrawText("Estrella (anti-horario)", 315, 25, 12, WHITE);
        DrawCircle(300, 55, 8, PURPLE);
        DrawText("Agujero (horario)", 315, 50, 12, WHITE);
        DrawCircle(300, 80, 8, GREEN);
        DrawText("Orbitas", 315, 75, 12, WHITE);

        if (cam.arrastrando) DrawCircle(GetMouseX(), GetMouseY(), 8, Color{ 255, 255, 255, 80 });

        EndDrawing();
    }

    CloseWindow();
    return 0;
}