#include <stdio.h>
#include <stdlib.h>
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
}

void mat_destroy(MAT* mat) { //Uvolnenie pamati (Звільнення пам'яті)
    free(mat->elem);
    free(mat);
}

int main() {
    char* f = "matica.bin";
    MAT* matica;
    unsigned int riadok, stlpec, var;
    printf("Rozmernost matice, format n n -> ");
    scanf("%i %i", &riadok, &stlpec);
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

    char subor_ano = mat_save(matica, f); //Vystup prvkov matice (Вивід елементів матриці)
    if (subor_ano) {
        printf("Matica je uspesne ulozena v subore.\n");
    }
    else {
        printf("Chyba pri ukladani matice do suboru.\n");
    }

    mat_destroy(matica); //Uvolnenie pamati (Звільнення пам'яті)

    return 0;
}
