#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ELEM(mat, i, j) ((mat)->elem[(i) * (mat)->stlpec + (j)])
#define FILE_OPEN_ERROR -1
#define HEADER_ERROR -2
#define READ_SIZE_ERROR -3
#define MATRIX_CREATE_ERROR -4
#define MATRIX_CSV_TABLE_ERROR -5

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
        return MATRIX_CREATE_ERROR;  // Chyba pri pridelovani pamate (Помилка виділення пам'яті)
    
    mat->riadok = riadok;
    mat->stlpec = stlpec;
    mat->elem = (float*)malloc(riadok * stlpec * sizeof(float));
    if (mat->elem == NULL) {
        free(mat);  // Chyba pri pridelovani pamate, uvolnenie pamate mat (Помилка виділення пам'яті, звільнення пам'яті mat)
        return MATRIX_CREATE_ERROR;
    }

    return mat;
}

MAT* mat_create_by_file(char* f) {
    FILE* file = fopen(f, "rb"); //Otvaranie suboru (Відкриття файлу)
    if (file == NULL) {
        return FILE_OPEN_ERROR;
    }

    // Kontrola prvych dvoch bajtov (Перевірка перших двох байтів)
    char header[2];
    if (fread(header, sizeof(char), 2, file) != 2 || header[0] != 'M' || header[1] != '1') {
        fclose(file);
        return HEADER_ERROR;
    }

    // Ziskanie poctu riadkov a stlpcov (Отримання кількості рядків і стовпців)
    unsigned int riadok, stlpec;
    if (fread(&riadok, sizeof(unsigned int), 1, file) != 1 || fread(&stlpec, sizeof(unsigned int), 1, file) != 1) {
        fclose(file);
        return READ_SIZE_ERROR;
    }

    MAT* matica = mat_create_with_type(riadok, stlpec); //(Ініціалізація матриці)
    if (matica == NULL) {
        fclose(file);
        return MATRIX_CREATE_ERROR;
    }

    // Citanie prvkov matice zo suboru (Зчитування елементів матриці з файлу)
    unsigned int num_elements = riadok * stlpec;
    if (fread(matica->elem, sizeof(float), num_elements, file) != num_elements) {
        fclose(file);
        free(matica->elem);
        mat_destroy(matica);
        return MATRIX_CREATE_ERROR;
    }

    fclose(file);
    return matica;
}

char mat_save(MAT* mat, char* f) {
    FILE* file = fopen(f, "wb");
    if (file == NULL) {
        return FILE_OPEN_ERROR; // Chyba pri otvarani suboru (Помилка при відкритті файлу)
    }

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

void mat_csv(const MAT* matica, const char* filename) { //https://matrixcalc.org/vectors.html
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        return MATRIX_CSV_TABLE_ERROR;
    }

    for (unsigned int i = 0; i < matica->riadok; i++) {
        for (unsigned int j = 0; j < matica->stlpec; j++) {
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
    for (unsigned int i = 0; i < mat->riadok; i++) {
        for (unsigned int j = 0; j < mat->stlpec; j++) {
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
    for (unsigned int i = 0; i < mat->riadok; i++) {
        for (unsigned int j = 0; j < mat->stlpec; j++) {
            float rand_cislo = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
            ELEM(mat, i, j) = rand_cislo;
        }
    }
}

void mat_print(MAT* mat) { //Vystup prvkov matice (Вивід елементів матриці)
    for (unsigned int i = 0; i < mat->riadok; i++) {
        for (unsigned int j = 0; j < mat->stlpec; j++) {
            printf("%6.2f\t", ELEM(mat, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

MAT* mat_multiply(MAT* mat1, MAT* mat2) { //Vypocet sucinu matic (Вираховування добутку матриць)
    unsigned int m = mat1->riadok;
    unsigned int n = mat1->stlpec;
    unsigned int p = mat2->stlpec;

    MAT* result = mat_create_with_type(m, p);
    if (result == NULL) {
        return MATRIX_CREATE_ERROR;
    }

    for (unsigned int i = 0; i < m; i++) { 
        for (unsigned int j = 0; j < p; j++) {
            float sum = 0.0;
            for (unsigned int k = 0; k < n; k++) {
                sum += ELEM(mat1, i, k) * ELEM(mat2, k, j);
            }
            ELEM(result, i, j) = sum;
        }
    }

    return result;
}

MAT* mat_subtract(MAT* mat1, MAT* mat2) { //Odpocitovanie matic
    if (mat1->riadok != mat2->riadok || mat1->stlpec != mat2->stlpec) {
        // Помилка: розміри матриць не співпадають
        return MATRIX_CREATE_ERROR;
    }

    MAT* result = mat_create_with_type(mat1->riadok, mat1->stlpec);
    if (result == NULL)
    {
        return MATRIX_CREATE_ERROR;
    }

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

    if (n == 1) {
        // Для матриці 1x1 детермінант дорівнює єдиному елементу
        return ELEM(mat, 0, 0);
    }

    float det = 0.0;
    int sign = 1;

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
    if (transposed == NULL) {
        return MATRIX_CREATE_ERROR;
    }

    for (unsigned int i = 0; i < mat->riadok; i++) {
        for (unsigned int j = 0; j < mat->stlpec; j++) {
            ELEM(transposed, j, i) = ELEM(mat, i, j);
        }
    }

    return transposed;
}

MAT* mat_orthogonalize(const MAT* mat) {
    unsigned int n = mat->riadok;

    MAT* orthogonalized = mat_create_with_type(mat->stlpec, mat->riadok);
    if (orthogonalized == NULL) {
        return MATRIX_CREATE_ERROR;
    }

    // Kopirovanie vstupnej matice do novej matice (Копіювання вхідної матриці до нової матриці)
    for (unsigned int i = 0; i < mat->riadok; i++) {
        for (unsigned int j = 0; j < mat->stlpec; j++) {
            ELEM(orthogonalized, i, j) = ELEM(mat, i, j);
        }
    }

    for (unsigned int i = 0; i < n; i++) {
        // Odcitanie projekcii vektorov na predchadzajuce bazove vektory
        for (unsigned int j = 0; j < i; j++) {
            float dot_product = 0.0;

            for (unsigned int k = 0; k < n; k++) {
                dot_product += ELEM(orthogonalized, k, i) * ELEM(orthogonalized, k, j);
            }

            for (unsigned int k = 0; k < n; k++) {
                ELEM(orthogonalized, k, i) -= dot_product * ELEM(orthogonalized, k, j);
            }
        }

        // Normalizacia vektorov (Нормалізація вектора)
        float norm = 0.0;
        for (unsigned int k = 0; k < n; k++) {
            norm += ELEM(orthogonalized, k, i) * ELEM(orthogonalized, k, i);
        }
        norm = sqrt(norm);

        if (norm != 0.0) {
            for (unsigned int k = 0; k < n; k++) {
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
    MAT* Q; //ortogonalna A
    MAT* R = mat_create_with_type(n, n); // Trojholnikova horna matica (Верхньо-трикутникова матриця)
    MAT* T = mat_create_with_type(n, n); // Transponovana Q

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
        T = mat_transpose(Q);
        A = mat_multiply(T, A);
        A = mat_multiply(A, Q);
        R = mat_multiply(T, A); //Vypocet iteracnych prvkov (Обчислення елементів ітерації)
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
    return 1;
}


int mat_error(int err)
{
    switch (err)
    {
    case -1:
        printf("Chyba: Otvorenie suboru zlyhalo\n");
        return 0;
    case -2:
        printf("Chyba: Chyba: Zapis typu matice v subore je nespravny\n");
        return 0;
    case -3:
        printf("Chyba: Nie je mozne urcit velkost matice zo suboru\n");
        return 0;
    case -4:
        printf("Chyba: Nepodarilo sa pridelit pamat pre maticu\n");
        return 0;
    case -5:
        printf("Chyba: Nepodarilo sa pridelit pamat pre zapis matice do tabulky\n");
        return 0;
    default:
        printf("Neznama chyba\n");
        return 0;
    }
}

int main() {
    srand(time(NULL));
    int i, chyba;
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
        if (matica == NULL) {
            printf("Chyba pri vytvarani matice.\n");
            return 1;
        }
        mat_random(matica);
        break;
    case 2:
        matica = mat_create_with_type(riadok, stlpec);
        if (matica == FILE_OPEN_ERROR || matica == HEADER_ERROR) {
            mat_error(matica);
            return 1;
        }
        mat_unit(matica);
        break;
    case 3:
        matica = mat_create_by_file(f, riadok);
        if (matica <= -3 && matica >= -8) {
            mat_error(matica);
            return 1;
        }
        if (matica == NULL) {
            printf("Chyba pri vytvarani matice zo suboru \n");
            return 1;
        }
        break;
    default: 
        printf("Nie je to spravna volba, skuste opat \n");
        mat_destroy(matica);
        return main();
    }

    float* coef = (float*)malloc((stlpec + 1) * sizeof(float));
    for (i = 0; i < stlpec + 1; i++)
    {
        coef[i] = 0.0;
    }

    mat_csv(matica, "matica.csv");

    // Vystup matice (Виведення матриці) 
    mat_print(matica);

    char subor_ano = mat_save(matica, f); 
    if (subor_ano) 
        printf("Matica je uspesne ulozena v subore.\n");
    
    else {
        printf("Chyba pri ukladani matice do suboru.\n");
        mat_destroy(matica);
        free(coef);
        return 0;
    }

    printf("Proces Ai+1 = Q*i.A.Q,\n");
    char polynom_ano = mat_characteristic_polynomial(matica, coef); 
    if (!polynom_ano) {
        printf("Chyba: Charakteristicky polynom sa nepodarilo vypocitat.\n");
        mat_destroy(matica);
        free(coef);
        return 0;
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
