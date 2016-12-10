
#include <iostream>
#include "mpi.h"
#include "limits.h"
#include "time.h"
#include "math.h"
#include "stdlib.h"
#include <vector>
#include <algorithm>

using namespace std;
using std::vector;

#define PART_A 0
#define PART_B 1
#define MANUAL -1
#define DEPTH 10

// Сортировка дефицитов по второму элементу
struct sort_def {
    bool operator()(const std::pair<int, int> &left, const std::pair<int, int> &right) {
        return left.second > right.second;
    }
};

class Graph {
    // Последняя итерация, на которой уменьшен разрез
    int last_improve;
    // Последняя пара обменянных вершин
    int last_a, last_b;
    // Размер графа
    int v_size, e_size;
    // Разбиения (текущее и лучшее)
    int *part;
    // Выгоды, массив времени обмена вершин
    int *when_moved, *gain;
    // Вектора дефицитов
    vector <pair<int, int> > def_a, def_b;
    // Нужно ли завершать итерацию
    bool is_finish;
public:
    // Матрица смежности графа
    int **edge;
    // Лучшее разбиение графа
    int *best_part;
    // Минимальный найденный разрез
    int min_cut;
    // Создание графа по размеру и вероятности ребра
    Graph(int size, int probty);
    Graph(int size, int *matr);
    // Удаление динамических структур графа
    ~Graph(void);
    // Печать данных для проверки корректности
    void print_matrix(void);
    void print_parts(void);
    void print_defs(void);
    void print_gains(void);
    void print_movtime(void);
    // Вычислеиие, пересчет и сортировка дефицитов
    void def_calc(void);
    void def_recalc(void);
    void def_balance(void);
    // Нужно ли завершать итерацию
    bool finish(void) { return is_finish; };
    void set_finish(bool finish) { is_finish = finish; }
    // Выбор частей случайным образом
    void random_parts(void);
    void choose_pair(int stage);
    // Шаг, на котором получены наиболее выгодные подмножества
    int max_step(void);
    // Их обмен на этом шаге
    void change_parts(int stage);
    // Разрез графа
    int cut(void);
    void copy_cut(int best_cut, int iter);
    // Вывод ответа
    void print_answer(void);
    // Сброс данных для параллельной версии
    void reset(void);
};
// Конструктор для рута
Graph::Graph(int size, int probty) {
    last_improve = 0;
    min_cut = INT_MAX;
    is_finish = false;
    unsigned i, j;
    int rem = size % 2;
    int half = (size + rem) / 2;
    e_size = 0;
    edge = new int*[size + rem];
    part = new int[size + rem];
    best_part = new int[size + rem];
    when_moved = new int[size + rem];
    gain = new int[half + 1];
    for (i = 0; i < half; i++) {
        part[i] = PART_A;
        part[i + half] = PART_B;
        gain[i] = 0;
    }
    gain[half] = 0;

    // Создаем матрицу
    for (i = 0; i < size + rem; i++) {
        edge[i] = new int[size + rem];
        for (j = 0; j < size + rem; j++)
            edge[i][j] = false;
    }

    // Генерируем матрицу
    if (probty != MANUAL) {
        int prob;
        for (i = 0; i < size; i++)
            for (j = i + 1; j < size; j++) {
                prob = rand() % 100;
                if (prob < probty) {
                    edge[i][j] = true;
                    edge[j][i] = true;
                    e_size++;
                }
            }
    }
    // Или заполняем вручную
    else {
        char buf;
        int value;
        for (i = 0; i < size; i++) {
            for (j = 0; j < size; j++) {
                cin >> buf;
                // Получаем числовое значение
                value = buf - '0';
                edge[i][j] = value;
                edge[j][i] = value;
            }
        }
    }
    // При необходимости добавляем фиктивную вершину
    if (rem == 1) {
        for (i = 0; i <= size; i++) {
            edge[size][i] = false;
            edge[i][size] = false;
        }
        // Теперь число вершин становится четным
        size++;
    }
    v_size = size;
}
// Конструктор для остальных процессов
Graph::Graph(int size, int *matr) {
    last_improve = 0;
    min_cut = INT_MAX;
    is_finish = false;
    unsigned i, j;
    int rem = size % 2;
    int half = (size + rem) / 2;
    e_size = 0;
    edge = new int*[size + rem];
    part = new int[size + rem];
    best_part = new int[size + rem];
    when_moved = new int[size + rem];
    gain = new int[half + 1];
    for (i = 0; i < half; i++) {
        part[i] = PART_A;
        part[i + half] = PART_B;
        gain[i] = 0;
    }
    gain[half] = 0;

    // Создаем матрицу
    for (i = 0; i < size + rem; i++) {
        edge[i] = new int[size + rem];
    }

    // Заполняем матрицу
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            edge[i][j] = matr[i*(size + rem) + j];
            if (edge[i][j]) e_size++;
        }
    }
    e_size /= 2;

    // При необходимости добавляем фиктивную вершину
    if (rem == 1) {
        for (i = 0; i <= size; i++) {
            edge[size][i] = false;
            edge[i][size] = false;
        }
        // Теперь число вершин становится четным
        size++;
    }
    v_size = size;
}

Graph::~Graph() {
    unsigned i;
    delete[] part;
    delete[] best_part;
    delete[] when_moved;
    delete[] gain;
    for (i = 0; i < v_size; i++)
        delete[] edge[i];
    delete edge;
}

// Печать матрицы смежности
void Graph::print_matrix() {
    unsigned i, j;
    cout << "There are " << e_size << " edges in the graph!" << endl;
    cout << "Here is the adjacency matrix:" << endl;
    for (i = 0; i < v_size; i++) {
        for (j = 0; j < v_size; j++)
            cout << edge[i][j];
        cout << endl;
    }
    cout << endl;
}

// Печать текущих частей разбиения
void Graph::print_parts() {
    unsigned i;
    cout << "The parts are:   ";
    for (i = 0; i < v_size; i++)
        cout << part[i];
    cout << endl;
}

// Печать пар дефицитов
void Graph::print_defs() {
    for (vector<pair<int, int> >::iterator it = def_a.begin(); it != def_a.end(); ++it)
        cout << it->second << " ";
    cout << endl;
    for (vector<pair<int, int> >::iterator it = def_b.begin(); it != def_b.end(); ++it)
        cout << it->second << " ";
    cout << endl;
}

// Печать массива выгод
void Graph::print_gains() {
    unsigned i;
    for (i = 0; i <= v_size / 2; i++)
        cout << gain[i] << " ";
    cout << endl;
}

void Graph::print_movtime() {
    unsigned i;
    for (i = 0; i < v_size; i++)
        cout << when_moved[i] << " ";
    cout << endl;
}

void Graph::print_answer() {
    unsigned i;
    cout << "The best partition is: ";
    for (i = 0; i < v_size; i++)
        cout << best_part[i];
    cout << endl;
    cout << "The last improvement has found on " << last_improve << " iteration!" << endl;
    cout << "The best cut contains " << min_cut << " vertices!" << endl;
}

void Graph::random_parts() {
    unsigned i;
    int dist = 0, v_random;
    for (i = 0; i < v_size; i++)
        part[i] = PART_A;
    while (dist < v_size / 2) {
        v_random = rand() % v_size;
        if (!part[v_random]) {
            part[v_random] = PART_B;
            dist++;
        }
    }
}
// Вычисление дефицитов
void Graph::def_calc() {
    int def;
    unsigned i, j;
    for (i = 0; i < v_size; i++) {
        when_moved[i] = 0;
    }

    for (i = 0; i < v_size; i++) {
        def = 0;
        // Вычисляем дефицит
        for (j = 0; j < v_size; j++)
            if (edge[i][j])
                if (part[i] == part[j]) def--; else def++;
        if (part[i] == PART_A)
            // И отправляем его в вектор, соответствующий вершине
            def_a.push_back(std::make_pair(i, def));
        else def_b.push_back(std::make_pair(i, def));
    }
}

// Пересчет дефицитов после обмена пары
void Graph::def_recalc() {
    // Проходим по вектору и пересчитываем дефицит каждой вершины
    for (vector<pair<int, int> >::iterator it = def_a.begin(); it != def_a.end(); ++it)
        it->second += 2 * edge[it->first][last_a] - 2 * edge[it->first][last_b];
    for (vector<pair<int, int> >::iterator it = def_b.begin(); it != def_b.end(); ++it)
        it->second += 2 * edge[it->first][last_b] - 2 * edge[it->first][last_a];
}

// Упорядочение дефицитов после обмена пары
void Graph::def_balance() {
    std::sort(def_a.begin(), def_a.end(), sort_def());
    std::sort(def_b.begin(), def_b.end(), sort_def());
}

// Обмен подмножеств
void Graph::change_parts(int stage) {
    unsigned i;
    if (stage > 0) {
        for (i = 0; i < v_size; i++)
            if (when_moved[i] <= stage)
                part[i] = !part[i];
    }
    else is_finish = true;
}

// Обмен пары вершин
void Graph::choose_pair(int stage) {
    unsigned i, j;
    int tgain, depth;
    int best_a, best_b;

    // Обычная оптимальная пара
    tgain = def_a[0].second + def_b[0].second - 2 * edge[def_a[0].first][def_b[0].first];
    best_a = 0;
    best_b = 0;

    // Граница на глубину поиска для экономии времени
    depth = min(DEPTH, v_size / 2 - stage + 1);

    // Соединенная пара может не являться оптимальной
    if (edge[def_a[0].first][def_b[0].first]) {
        for (i = 0; i < depth; i++)
            for (j = 0; j < depth; j++)
                // Если нашли пару лучше уже найденной
                if (def_a[i].second + def_b[j].second - 2 * edge[def_a[i].first][def_b[j].first] > tgain) {
                    best_a = i;
                    best_b = j;
                    tgain = def_a[i].second + def_b[j].second - 2 * edge[def_a[i].first][def_b[j].first];
                }
    }
    // Оптимальная пара
    last_a = def_a[best_a].first;
    last_b = def_b[best_b].first;
    // Удаляем из векторов пару вершин, которые здесь обмениваем
    def_a.erase(def_a.begin() + best_a);
    def_b.erase(def_b.begin() + best_b);
    // И указываем данные об обмене
    when_moved[last_a] = stage;
    when_moved[last_b] = stage;
    gain[stage] = gain[stage - 1] + tgain;
}

// Вычисление шага, на котором выгода от обмена максимальна
int Graph::max_step() {
    unsigned i;
    int gmax = 0, stage = 0;
    for (i = 1; i <= v_size / 2; i++) {
        if (gain[i] > gmax) {
            stage = i;
            gmax = gain[i];
        }
    }
    return stage;
}

// Возвращает величину разреза
int Graph::cut() {
    unsigned i, j;
    unsigned cut = 0;
    for (i = 0; i < v_size; i++)
        for (j = i + 1; j < v_size; j++)
            if ((edge[i][j]) && (part[i] != part[j]))
                cut++;
    return cut;
}

// Если текущий разрез лучше минимального
void Graph::copy_cut(int best_cut, int iter) {
    unsigned i;
    last_improve = iter;
    min_cut = best_cut;
    for (i = 0; i < v_size; i++)
        best_part[i] = part[i];
}

void Graph::reset() {
    int half = v_size / 2, i;
    for (i = 0; i < half; i++) {
        part[i] = PART_A;
        part[i + half] = PART_B;
        gain[i] = 0;
    }
    gain[half] = 0;
    last_improve = 0;
    min_cut = INT_MAX;
    is_finish = false;
}

// Сколько итераций выполнять каждому процессу
int iter_per_process(int rank, int total_iter) {
    if (rank > total_iter) return 1;
    if (total_iter % rank == 0) return total_iter / rank;
    else return (total_iter / rank + 1);
}

int main(int argc, char** argv) {
    int process_min, best_result, best_process = 0;
    int v_size, e_size, prob, probty, rem, t_cut;
    int rank, my_rank;
    unsigned i, j, iter, iters;
    int *final_part;
    int *matrix;

    double start, finish, stime, ptime, accel;
    const int root = 0;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    srand(time(NULL)*time(NULL)*my_rank*my_rank*my_rank);
    // Первый аргумент — число итераций, второй — число вершин
    // третий - вероятность (в процентах) появления каждого ребра
    if (my_rank == root) {
        // Значения по умолчанию
        iter = 50;
        probty = 50;
        v_size = 6;
        // v_size = atoi(argv[1]);
        //iter = atoi(argv[2]);
        //probty = atoi(argv[3]);
        rem = v_size % 2;
        start = MPI_Wtime();
        Graph G(v_size, probty);
        for (j = 1; j <= iter; j++) {
            G.set_finish(false);
            G.random_parts();
            // Пока можно уменьшить разрез
            while (!G.finish()) {
                G.def_calc();
                G.def_balance();
                // Печать нужна для проверки корректности алгоритма
                for (i = 1; i <= (v_size + 1) / 2; i++) {
                    G.choose_pair(i);
                    G.def_recalc();
                    G.def_balance();
                }
                G.change_parts(G.max_step());
            }
            t_cut = G.cut();
            if (t_cut < G.min_cut) { G.copy_cut(t_cut, j); }
        }
        finish = MPI_Wtime();
        G.print_matrix();
        cout << "The best partition is: ";
        for (i = 0; i < v_size + 1; i++)
            cout << G.best_part[i];
        cout << endl;
        cout << "The best cut contains " << G.min_cut << " edges!" << endl;
        stime = finish - start;
        printf("Serial time: %f\n", stime);
        // Удаляем данные последовательной версии
        G.reset();
        // Параллельная версия
        if (rank != 1) {
            start = MPI_Wtime();
            MPI_Barrier(MPI_COMM_WORLD);
            iters = iter_per_process(rank, iter);
            MPI_Bcast(&v_size, 1, MPI_INT, root, MPI_COMM_WORLD);
            MPI_Bcast(&iters, 1, MPI_INT, root, MPI_COMM_WORLD);
            rem = v_size % 2;
            matrix = new int[(v_size + rem)*(v_size + rem)];
            final_part = new int[v_size + rem];
            for (i = 0; i < v_size + rem; i++) {
                for (j = 0; j < v_size + rem; j++)
                    matrix[i*(v_size + rem) + j] = G.edge[i][j];
            }

            // Отправляем матрицу другим процессам
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(matrix, (v_size + rem)*(v_size + rem), MPI_INT, root, MPI_COMM_WORLD);
            // Алгоритм разделения графа
            for (j = 1; j <= iters; j++) {
                G.set_finish(false);
                G.random_parts();
                // Пока можно уменьшить разрез
                while (!G.finish()) {
                    G.def_calc();
                    G.def_balance();
                    // Печать нужна для проверки корректности алгоритма
                    for (i = 1; i <= (v_size + 1) / 2; i++) {
                        G.choose_pair(i);
                        G.def_recalc();
                        G.def_balance();
                    }
                    G.change_parts(G.max_step());
                }
                t_cut = G.cut();
                if (t_cut < G.min_cut) { G.copy_cut(t_cut, j); }
            }

            for (i = 0; i < v_size + rem; i++)
                final_part[i] = G.best_part[i];
            best_result = G.min_cut;
            for (i = 1; i < rank; i++) {
                MPI_Recv(&process_min, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
                if (process_min < best_result) {
                    best_process = i;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(&best_process, 1, MPI_INT, root, MPI_COMM_WORLD);
            if (best_process != root) {
                MPI_Recv(final_part, v_size + rem, MPI_INT, best_process, 1, MPI_COMM_WORLD, &status);
            }
           /* cout << "The best partition is: ";
            for (i = 0; i < v_size + rem; i++)
                cout << final_part[i];
            cout << endl;
            cout << "The best cut contains " << best_result << " edges!" << endl;*/
            finish = MPI_Wtime();
            ptime = finish - start;
            printf("Parallel time: %f\n", ptime);
            accel = stime / ptime;
            printf("The acceleration is: %f\n", accel);
            delete[] final_part;
            delete[] matrix;
        }
    }
    else {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&v_size, 1, MPI_INT, root, MPI_COMM_WORLD);
        rem = v_size % 2;
        MPI_Bcast(&iters, 1, MPI_INT, root, MPI_COMM_WORLD);
        matrix = new int[(v_size + rem)*(v_size + rem)];
        final_part = new int[v_size + rem];
        // Принимаем матрицу смежности
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(matrix, (v_size + rem)*(v_size + rem), MPI_INT, root, MPI_COMM_WORLD);
        Graph G = Graph(v_size, matrix);
        // Алгоритм разделения графа
        for (j = 1; j <= iters; j++) {
            G.set_finish(false);
            G.random_parts();
            // Пока можно уменьшить разрез
            while (!G.finish()) {
                G.def_calc();
                G.def_balance();
                // Печать нужна для проверки корректности алгоритма
                for (i = 1; i <= (v_size + 1) / 2; i++) {
                    G.choose_pair(i);
                    G.def_recalc();
                    G.def_balance();
                }
                G.change_parts(G.max_step());
            }
            t_cut = G.cut();
            if (t_cut < G.min_cut) { G.copy_cut(t_cut, j); }
        }
        process_min = G.min_cut;
        MPI_Send(&process_min, 1, MPI_INT, root, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&best_process, 1, MPI_INT, root, MPI_COMM_WORLD);
        if (best_process == my_rank) {
            MPI_Send(G.best_part, v_size + rem, MPI_INT, root, 1, MPI_COMM_WORLD);
        }
        delete[] final_part;
        delete[] matrix;
    }
    MPI_Finalize();
    return 0;
}