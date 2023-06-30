#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ELEM(mat, i, j) ((mat)->elem[(i) * (mat)->stlpec + (j)])
#define POLYNOM_ANO 0
#define POLYNOM_NIE 1
#define FILE_OPEN_ERROR -1
#define MATRIX_CSV_TABLE_ERROR -2

//Ak v komentaroch vidite nejasne znaky, je to ukrajincina
//https://matrixcalc.org/vectors.html
typedef struct {
    unsigned int riadok;
    unsigned int stlpec;
    float* elem;
} MAT;

void mat_destroy(MAT* mat) { //Uvolnenie pamati (Звільнення пам'яті)
    free(mat->elem);
    free(mat);
}

MAT* mat_create_with_type(unsigned int riadok, unsigned int stlpec) {
    MAT* mat = (MAT*)malloc(sizeof(MAT));
    if (mat == NULL) 
        return NULL;  // Chyba pri pridelovani pamate (Помилка виділення пам'яті)
    
    mat->riadok = riadok;
    mat->stlpec = stlpec;
    mat->elem = (float*)malloc(riadok * stlpec * sizeof(float));
    if (mat->elem == NULL) {
        free(mat);  // Chyba pri pridelovani pamate, uvolnenie pamate mat (Помилка виділення пам'яті, звільнення пам'яті mat)
        return NULL;
    }

    return mat;
}

MAT* mat_create_by_file(char* f) {
    unsigned int riadok, stlpec, num_elements;
    FILE* file = fopen(f, "rb"); //Otvaranie suboru (Відкриття файлу)
    if (file == NULL) 
        return NULL;

    // Kontrola prvych dvoch bajtov (Перевірка перших двох байтів)
    char header[2];
    if (fread(header, sizeof(char), 2, file) != 2 || header[0] != 'M' || header[1] != '1') {
        fclose(file);
        return NULL;
    }

    // Ziskanie poctu riadkov a stlpcov (Отримання кількості рядків і стовпців)
    if (fread(&riadok, sizeof(unsigned int), 1, file) != 1 || fread(&stlpec, sizeof(unsigned int), 1, file) != 1)
    {
        fclose(file);
        return NULL;
    }

    MAT* matica = mat_create_with_type(riadok, stlpec); //(Ініціалізація матриці)
    if (matica == NULL) 
    {
        fclose(file);
        return NULL;
    }

    // Citanie prvkov matice zo suboru (Зчитування елементів матриці з файлу)
    num_elements = riadok * stlpec;
    if (fread(matica->elem, sizeof(float), num_elements, file) != num_elements) 
    {
        fclose(file);
        free(matica->elem);
        mat_destroy(matica);
        return NULL;
    }

    fclose(file);
    return matica;
}

char mat_save(MAT* mat, char* f) {
    FILE* file = fopen(f, "wb");
    if (file == NULL) 
        return FILE_OPEN_ERROR; // Chyba pri otvarani suboru (Помилка при відкритті файлу)

    // Zaznamenanie nazvu(header)
    char header[] = { 'M', '1' };
    fwrite(header, sizeof(char), 2, file);

    // Zaznamenanie poctu riadkov a stlpcov (Запис к-сті радків та стовпців)
    fwrite(&(mat->riadok), sizeof(unsigned int), 1, file);
    fwrite(&(mat->stlpec), sizeof(unsigned int), 1, file);

    // Zaznamenavanie prvkov matice (Запис елементів матриці)
    fwrite(mat->elem, sizeof(float), mat->riadok * mat->stlpec, file);

    fclose(file);
    return 1;
}

char mat_csv(const MAT* matica, const char* filename) { //https://matrixcalc.org/vectors.html
    unsigned int i, j;
    
    FILE* file = fopen(filename, "w");
    if (file == NULL) 
        return MATRIX_CSV_TABLE_ERROR;

    for (i = 0; i < matica->riadok; i++) {
        for (j = 0; j < matica->stlpec; j++) {
            fprintf(file, "%.2f", ELEM(matica, i, j));
            if (j < matica->stlpec - 1) {
                fprintf(file, ",");
            }
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

void mat_unit(MAT* mat) { //Generovanie jednotkovej matice (Генерування одиничної матриці)
    unsigned int i, j;

    for (i = 0; i < mat->riadok; i++) {
        for (j = 0; j < mat->stlpec; j++) {
            if (i == j) {
                ELEM(mat, i, j) = 1.0;
            }
            else {
                ELEM(mat, i, j) = 0.0;
            }
        }
    }
}

void mat_random(MAT* mat) { //Nahodne generovanie prvkov matice (Випадкова генерація елементів матриці)
    unsigned int i, j;
    float rand_cislo;

    for (i = 0; i < mat->riadok; i++) {
        for (j = 0; j < mat->stlpec; j++) {
            rand_cislo = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
            ELEM(mat, i, j) = rand_cislo;
        }
    }
}

void mat_print(MAT* mat) { //Vystup prvkov matice (Вивід елементів матриці)
    unsigned int i, j;
    for (i = 0; i < mat->riadok; i++) {
        for (j = 0; j < mat->stlpec; j++) {
            printf("%6.2f\t", ELEM(mat, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

MAT* mat_multiply(MAT* mat1, MAT* mat2) { //Vypocet sucinu matic (Вираховування добутку матриць)
    unsigned int i, j, k;
    unsigned int m = mat1->riadok;
    unsigned int n = mat1->stlpec;
    unsigned int p = mat2->stlpec;

    MAT* result = mat_create_with_type(m, p);
    if (result == NULL) 
        return NULL;
    
    for (i = 0; i < m; i++) { 
        for (j = 0; j < p; j++) {
            float sum = 0.0;
            for (k = 0; k < n; k++) {
                sum += ELEM(mat1, i, k) * ELEM(mat2, k, j);
            }
            ELEM(result, i, j) = sum;
        }
    }

    return result;
}

MAT* mat_subtract(MAT* mat1, MAT* mat2) { //Odpocitovanie matic
    if (mat1->riadok != mat2->riadok || mat1->stlpec != mat2->stlpec) 
        return NULL;
    

    MAT* result = mat_create_with_type(mat1->riadok, mat1->stlpec);
    if (result == NULL)
        return NULL;

    for (unsigned int i = 0; i < mat1->riadok; i++) {
        for (unsigned int j = 0; j < mat1->stlpec; j++) {
            ELEM(result, i, j) = ELEM(mat1, i, j) - ELEM(mat2, i, j);
        }
    }

    return result;
}

float mat_determinant(MAT* mat) {
    unsigned int n = mat->riadok;
    int i, j, k;
    float det = 0.0;
    int sign = 1;

    if (n == 1)
        return ELEM(mat, 0, 0);

    for (j = 0; j < n; j++) {
        MAT* submatrix = mat_create_with_type(n - 1, n - 1);

        for (i = 1; i < n; i++) {
            unsigned int subi = 0;
            for (k = 0; k < n; k++) {
                if (k != j) {
                    ELEM(submatrix, subi, i - 1) = ELEM(mat, i, k);
                    subi++;
                }
            }
        }

        float subdet = mat_determinant(submatrix);
        det += sign * ELEM(mat, 0, j) * subdet;

        sign = -sign;

        mat_destroy(submatrix);
    }

    return det;
}

MAT* mat_transpose(const MAT* mat) { //Transponovanie matice (Транспонування матриці)
    MAT* transposed = mat_create_with_type(mat->stlpec, mat->riadok);
    unsigned int i, j;

    if (transposed == NULL) 
        return NULL;
    

    for (i = 0; i < mat->riadok; i++) {
        for (j = 0; j < mat->stlpec; j++) {
            ELEM(transposed, j, i) = ELEM(mat, i, j);
        }
    }

    return transposed;
}

MAT* mat_orthogonalize(const MAT* mat) {
    unsigned int n = mat->riadok, i, j, k;

    MAT* orthogonalized = mat_create_with_type(mat->stlpec, mat->riadok);
    if (orthogonalized == NULL) 
        return NULL;

    // Kopirovanie vstupnej matice do novej matice (Копіювання вхідної матриці до нової матриці)
    for (i = 0; i < mat->riadok; i++) {
        for (j = 0; j < mat->stlpec; j++) {
            ELEM(orthogonalized, i, j) = ELEM(mat, i, j);
        }
    }

    for (i = 0; i < n; i++) {
        // Odcitanie projekcii vektorov na predchadzajuce bazove vektory
        for (j = 0; j < i; j++) {
            float dot_product = 0.0;

            for (k = 0; k < n; k++) {
                dot_product += ELEM(orthogonalized, k, i) * ELEM(orthogonalized, k, j);
            }

            for (k = 0; k < n; k++) {
                ELEM(orthogonalized, k, i) -= dot_product * ELEM(orthogonalized, k, j);
            }
        }

        // Normalizacia vektorov (Нормалізація вектора)
        float norm = 0.0;
        for (k = 0; k < n; k++) {
            norm += ELEM(orthogonalized, k, i) * ELEM(orthogonalized, k, i);
        }
        norm = sqrt(norm);

        if (norm != 0.0) {
            for (k = 0; k < n; k++) {
                ELEM(orthogonalized, k, i) /= norm;
            }
        }
    }

    return orthogonalized;
}

/*Hladanie charakteristickeho polynomu sa realizuje prostrednictvom principu "Schurovho rozkladu", kde vzorec je Ai+1 = Q*i.A.Q,
kde Ai+1 je dalsia matica, ktora je podobnejsia hornej trojuholnikovej matici ako Ai, ale tolko do pevneho i, Q je ortogonalizovana matica A a Q* je
transponovana matica Q. Ide o to, ze koncova A bude mat vlastne cisla originalnej A. 
Iteracie Ai+1 = Q*i.A.Q su potrebne na to, aby sa vlastne hodnoty diagonaly A co najviac priblizili realnym vlastnym hodnotam A, ak vascsie
A je podobne do trojholnikovej, tak prvky na jej diagonale je blizsie do vlastnych hodnot zaciatocnej A*/
char mat_characteristic_polynomial(MAT* mat, float* coef) {
    int i, j, ii, jj;
    unsigned int n = mat->riadok;
    MAT* A = mat_create_with_type(n, n); // Kopia vstupnej matice (Копія вхідної матриці)
    if (A == NULL)
        return POLYNOM_NIE;
    MAT* Q; //ortogonalna A
    MAT* R = mat_create_with_type(n, n); // Trojholnikova horna matica (Верхньо-трикутникова матриця)
    if (R == NULL)
        return POLYNOM_NIE;
    MAT* T = mat_create_with_type(n, n); // Transponovana Q
    if (T == NULL)
        return POLYNOM_NIE;

    // Kopirovanie vstupnej matice do novej matice (Копіювання вхідної матриці до нової матриці)
    for (unsigned int i = 0; i < mat->riadok; i++) {
        for (unsigned int j = 0; j < mat->stlpec; j++) {
            ELEM(A, i, j) = ELEM(mat, i, j);
        }
    }

    //Vypocet poctu nul v hornej trojuholnikovej diagonale (pocet prvkov pod diagonalou) (Обчислення кількості нулів у верхньо трикутній діагоналі(к-сть елементів під діагоналлю))
    int k = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i > j) {
                k++;
            }
        }
    }
    //printf("\n k = %i\n", k);
    float presnost = 0.01; // Prahova hodnota presnosti (Порогове значення точності)

    // Vykonavanie procesu dosahovania presnosti (Виконання процесу досягнення точності)
    for (i = 0; ; i++)
    {

        Q = mat_orthogonalize(A);
        if (Q == NULL)
            return POLYNOM_NIE;
        T = mat_transpose(Q);
        if (T == NULL)
            return POLYNOM_NIE;
        A = mat_multiply(T, A);
        if (A == NULL)
            return POLYNOM_NIE;
        A = mat_multiply(A, Q);
        if (A == NULL)
            return POLYNOM_NIE;
        R = mat_multiply(T, A); //Vypocet iteracnych prvkov (Обчислення елементів ітерації)
        if (R == NULL)
            return POLYNOM_NIE;
        int je_podobne = 0;

        for (ii = 0; ii < n; ii++) { //pouzivam tu R na overenie presnosti, pretoze zo A cyklus stava velmi velky(prakticke nekonecny)
            for (jj = 0; jj < n; jj++) {
                if (ii > jj) {
                    if (fabs(ELEM(R, ii, jj)) < presnost)
                        je_podobne++; //Sucet poctu prvkov R pod diagonalou, ktore su mensie ako presnost (Кількість елементів матриці R рід діагоналлю, які менші за коефіцієнт точності)
                }
            }
        }
        if (je_podobne == k) 
            break; // Zastavenie cyklu po dosiahnuti pozadovanej presnosti (Зупиняємо цикл, якщо досягнута потрібна точність)
        
    }
    printf("Iteracia %i je najlepsia, pri ktorej matica R je najviac podobna hornej trojuholnikovej matici s presnostou \"nulovych\" prvkov do %f absolutnej hodnoty,\nA\n", i, presnost);
    mat_print(A);
    printf("R\n"); //zaujimave, ze na diagonale R je vlaste cisla A, ale oni su vsetke kladne!!!
    mat_print(R);

    coef[n] = 1.0;

    for (int i = 0; i < n; i++) { //Metoda podla Viethovho vzorca
        for (int j = n - i - 1; j < n; j++) {
            coef[j] = coef[j] - ELEM(A, i, i) * coef[j + 1];
        }
    }

    mat_destroy(A);
    mat_destroy(T);
    mat_destroy(Q);
    mat_destroy(R);
    return POLYNOM_ANO;
}

int main() {
    srand(time(NULL));
    int i;
    char polynom_ano, subor_ano;
    char* f = "matica.bin";
    MAT* matica;
    unsigned int riadok, stlpec, var;

    printf("Rozmernost matice, matica je kvadraticka -> ");
    scanf("%i", &riadok);
    stlpec = riadok;
    if (riadok != stlpec)
    {
        printf("Matica nie je kvadraticka!\n");
        return main();
    }

    printf("1 - rand generacia matice, 2 - jednotkova diagonalna matica, 3 - nacitanie zo suboru -> ");
    scanf("%i", &var);

    switch (var) { //moznosti generovania matic (варіанти генерування матриці)
    case 1:
        matica = mat_create_with_type(riadok, stlpec);
        if (matica == NULL)
            return -1;
        mat_random(matica);
        break;
    case 2:
        matica = mat_create_with_type(riadok, stlpec);
        if (matica == NULL)
            return -1;
        mat_unit(matica);
        break;
    case 3:
        matica = mat_create_by_file(f);
        if (matica == NULL)
            return -1;
        break;
    default: 
        printf("Nie je to spravna volba, skuste opat \n");
        mat_destroy(matica);
        return main();
    }

    // Vystup matice (Виведення матриці) 
    mat_print(matica);

    float* coef = (float*)malloc((stlpec + 1) * sizeof(float));
    for (i = 0; i < stlpec + 1; i++)
    {
        coef[i] = 0.0;
    }

    subor_ano = mat_csv(matica, "matica.csv");
    if (subor_ano == MATRIX_CSV_TABLE_ERROR)
    {
        printf("Chyba pri ukladani matice do tabulky.\n");
        mat_destroy(matica);
        free(coef);
        return -1;
    }

    subor_ano = mat_save(matica, f); 
    if (subor_ano == FILE_OPEN_ERROR)
    {
        printf("Chyba pri ukladani matice do suboru.\n");
        mat_destroy(matica);
        free(coef);
        return -1;
    }

    printf("Matica je uspesne ulozena v subore.\n");

    printf("Proces Ai+1 = Q*i.A.Q,\n");
    polynom_ano = mat_characteristic_polynomial(matica, coef); 
    if (polynom_ano == POLYNOM_NIE) {
        printf("Chyba: Charakteristicky polynom sa nepodarilo vypocitat.\n");
        mat_destroy(matica);
        free(coef);
        return -1;
    }

    printf("Polynom: ");
    for (i = 0; i < stlpec+1; i++)
    {
        printf("%6.5f ", coef[i]);
    }
    printf("\n");

    free(coef);
    mat_destroy(matica); //Uvolnenie pamati (Звільнення пам'яті)

    return 0;
}
