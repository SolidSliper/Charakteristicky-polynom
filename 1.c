#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ELEM(mat, i, j) ((mat)->elem[(i) * (mat)->stlpec + (j)])

//Ak v komentaroch vidite nejasne znaky, je to ukrajincina
typedef struct {
    unsigned int riadok;
    unsigned int stlpec;
    float* elem;
} MAT;



MAT* mat_create_with_type(unsigned int riadok, unsigned int stlpec) {
    MAT* mat = (MAT*)malloc(sizeof(MAT));
    if (mat == NULL) {
        return NULL;  // Chyba pri pridelovani pamate (Помилка виділення пам'яті)
    }

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
    FILE* file = fopen(f, "rb"); //Otvaranie suboru (Відкриття файлу)
    if (file == NULL) {
        printf("Chyba pri otvarani suboru.\n");
        return NULL;
    }

    // Kontrola prvych dvoch bajtov (Перевірка перших двох байтів)
    char header[2];
    if (fread(header, sizeof(char), 2, file) != 2 || header[0] != 'M' || header[1] != '1') {
        printf("Nespravny format suboru.\n");
        fclose(file);
        return NULL;
    }

    // Ziskanie poctu riadkov a stlpcov (Отримання кількості рядків і стовпців)
    unsigned int riadok, stlpec;
    if (fread(&riadok, sizeof(unsigned int), 1, file) != 1 || fread(&stlpec, sizeof(unsigned int), 1, file) != 1) {
        printf("Chyba pri citani rozmerov matice.\n");
        fclose(file);
        return NULL;
    }

    MAT* matica = mat_create_with_type(riadok, stlpec); //(Ініціалізація матриці)
    if (matica == NULL) {
        printf("Chyba pri vytvarani matice.\n");
        fclose(file);
        return NULL;
    }

    // Citanie prvkov matice zo suboru (Зчитування елементів матриці з файлу)
    unsigned int num_elements = riadok * stlpec;
    if (fread(matica->elem, sizeof(float), num_elements, file) != num_elements) {
        printf("Chyba pri citani prvkov matice. \n");
        fclose(file);
        free(matica->elem);
        free(matica);
        return NULL;
    }

    fclose(file);
    return matica;
}



char mat_save(MAT* mat, char* f) {
    FILE* file = fopen(f, "wb");
    if (file == NULL) {
        return 0; // Chyba pri otvarani suboru (Помилка при відкритті файлу)
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
    srand(time(NULL));

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



void mat_destroy(MAT* mat) { //Uvolnenie pamati (Звільнення пам'яті)
    free(mat->elem);
    free(mat);
}



MAT* mat_multiply(MAT* mat1, MAT* mat2) { //Vypocet sucinu matic (Вираховування добутку матриць)
    unsigned int m = mat1->riadok;
    unsigned int n = mat1->stlpec;
    unsigned int p = mat2->stlpec;

    MAT* result = mat_create_with_type(m, p);
    if (result == NULL) {
        printf("Chyba: Nepodarilo sa vytvorit maticu pre vysledok.\n");
        return NULL;
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



MAT* mat_transpose(const MAT* mat) { //Transponovanie matice (Транспонування матриці)
    MAT* transposed = mat_create_with_type(mat->stlpec, mat->riadok);
    if (transposed == NULL) {
        printf("Chyba pri vytvarani transponovanej matice.\n");
        return NULL;
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
        printf("Chyba pri vytvarani ortogonalizovanej matice.\n");
        return NULL;
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
transponovana matica Q. Ide o to, ze mame aj vzorec R = Q*.A, kde R je horna trojuholnikova matica s vlastnymi hodnotami A na diagonale. 
Iteracie Ai+1 = Q*i.A.Q su potrebne na to, aby sa vlastne hodnoty diagonaly R co najviac priblizili realnym vlastnym hodnotam A*/
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
    printf("\n k = %i\n", k);
    float presnost = 0.01; // Prahova hodnota presnosti (Порогове значення точності)

    // Vykonavanie procesu dosahovania presnosti (Виконання процесу досягнення точності)
    for (i = 0; ; i++)
    {

        Q = mat_orthogonalize(A);
        T = mat_transpose(Q);
        A = mat_multiply(T, A);
        A = mat_multiply(A, Q);
        R = mat_multiply(T, A); //Vypocet iteracnych prvkov (Обчислення елементів ітерації)
        //mat_print(R);
        int je_podobne = 0;

        for (ii = 0; ii < n; ii++) {
            for (jj = 0; jj < n; jj++) {
                if (ii > jj) {
                    if (fabs(ELEM(R, ii, jj)) < presnost)
                        je_podobne++; //Sucet poctu prvkov R pod diagonalou, ktore su mensie ako presnost (Кількість елементів матриці R рід діагоналлю, які менші за коефіцієнт точності)
                }
                //coef[ii] = ELEM(R, ii, ii);
            }
        }
        if (je_podobne == k) {
            break; // Zastavenie cyklu po dosiahnuti pozadovanej presnosti (Зупиняємо цикл, якщо досягнута потрібна точність)
        }
    }
    R = mat_multiply(T, A);
    printf("Iteracia %i je najlepsia, pri ktorej matica R je najviac podobna hornej trojuholnikovej matici s presnostou \"nulovych\" prvkov do %f absolutnej hodnoty,\nA\n", i, presnost);
    mat_print(A);
    for (i = 0; i < n; i++) // Zapis charakterestickeho polynomu
    {
        coef[i] = ELEM(R, i, i);
    }
    printf("R\n");
    mat_print(R);

    mat_destroy(A);
    mat_destroy(T);
    mat_destroy(Q);
    mat_destroy(R);
    return 1;
}

int main() {
    int i;
    char* f = "matica.bin";
    MAT* matica;
    unsigned int riadok, stlpec, var;
    printf("Rozmernost matice, matica je kvadraticka -> ");
    scanf("%i", &riadok);
    stlpec = riadok;
    float* coef = (float*)malloc((stlpec + 1) * sizeof(float));
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
        if (matica == NULL) {
            printf("Chyba pri vytvarani matice \n");
            return 1;
        }
        mat_unit(matica);
        break;
    case 3:
        matica = mat_create_by_file(f);
        if (matica == NULL) {
            printf("Chyba pri vytvarani matice zo suboru \n");
            return 1;
        }
        break;
    default: 
        printf("Nie je to spravna volba, skuste opat \n");
        return main();
    }

    // Vystup matice (Виведення матриці) 
    mat_print(matica);

    char subor_ano = mat_save(matica, f); 
    if (subor_ano) {
        printf("Matica je uspesne ulozena v subore.\n");
    }
    else {
        printf("Chyba pri ukladani matice do suboru.\n");
        return 0;
    }

    printf("Proces Ai+1 = Q*i.A.Q,\n");
    mat_characteristic_polynomial(matica, coef);

    char polynom_ano = mat_characteristic_polynomial(matica, coef); 
    if (!polynom_ano) {
        printf("Chyba: Charakteristicky polynom sa nepodarilo vypocitat.\n");
        return 0;
    }
    printf("Polynom: ");
    for (i = 0; i < stlpec; i++)
    {
        printf("%6.5f ", coef[i]);
    }

    free(coef);
    mat_destroy(matica); //Uvolnenie pamati (Звільнення пам'яті)

    return 0;
}
