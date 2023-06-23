# Charakteristicky polynom

Funkcia funguje dobre, nájde najbližšie možné hodnoty vlastných čísel za predpokladu, že tieto čísla sú reálne, 
pri imaginárnych číslach jednoducho vráti určité reálne číslo, rovnaký princíp hľadania QR alebo Q*AQ

Ak chcete skontrolovať vlastné čísla, použite webovú stránku https://matrixcalc.org/vectors.html, stačí skopírovať maticu z tabuľky matica.csv

Uloha 9.1.14: Okrem ˇstandardnej funkcionality naprogramujte aj funkciu
• char mat_characteristic_pollynomial(MAT *mat, float* coef);
ktor´a vypoˇc´ıta koeficienty charakteristick´eho polyn´omu matice mat. Cez argu-
ment coef funkcia prevezme adresu pol’a d´lˇzky n + 1 kde n je typ matice mat,
do ktor´eho uloˇz´ı zisten´e koeficienty charakteristick´eho polyn´omu. V pr´ıpade,
ˇze funkcia z nejak´eho fundament´alneho dˆovodu nevie charakteristick´y polyn´om
vypoˇc´ıtat’, st’aˇzuje sa cez n´avratov´u hodnotu
