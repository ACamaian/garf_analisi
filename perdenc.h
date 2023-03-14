#ifdef __cplusplus
extern "C" {
#endif

typedef struct Matrl {
	char nome[32];
	int ne; //numero elementi (MAX 5)
	int e[5][3]; //elementi (numero atomico + numero di massa (-1 per A medio da distribuzione naturale) + numero stechiometrico)
	double r; //densità in mg/cm². -1 se è un gas composto o un elemento gassoso, 0 se è un elemento solido. Valore se è un composto solido.
	double P; //pressione in mbar (SOLO GAS, viene ignorata per i solidi);
	double T; //temperatura in K (SOLO GAS, viene ignorata per i solidi);
} MATE;

//restituisce l'energia dopo l'assorbitore a partire dall'energia prima
double eloss(double ein /*in MeV*/,int Zp,int Ap,MATE *mat,double thickness /*in µm*/,int form);

//restituisce l'energia prima dell'assorbitore a partire dall'energia dopo
double eorig(double epass /*in MeV*/,int Zp,int Ap,MATE *mat,double thickness /*in µm*/,int form);

//restituisce il range data l'energia iniziale
double e2range(double ein /*in MeV*/,int Zp,int Ap,MATE *mat,int form);

//restituisce l'energia iniziale a partire dal range
double range2e(double range /*in mg/cm²*/,int Zp,int Ap,MATE *mat,int form);

//restituisce lo spessore dell'assorbitore in mg/cm² date l'energia iniziale e l'energia finale
double e2spes(double ein /*in MeV*/,int Zp,int Ap,MATE *mat,double epass /*in MeV*/,int form);

//restituisce l'energia iniziale data l'energia rilasciata e lo spessore
double spes2e(double el /*in MeV*/,int Zp,int Ap,MATE *mat,double thickness /*in µm*/,int form);

#ifdef __cplusplus
}
#endif
