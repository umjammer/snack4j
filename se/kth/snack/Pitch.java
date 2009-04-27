/* 
 * Copyright (C) 2000-2004 Kåre Sjölander <kare@speech.kth.se>
 * Copyright (C) 1997 Philippe Langlais <felipe@speech.kth.se>
 *
 * This file is part of the Snack Sound Toolkit.
 * The latest version can be found at http://www.speech.kth.se/snack/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

package se.kth.snack;

public class Pitch implements Command {

    /* 
     * LES PARAMETRES GLOBAUX DU DETECTEUR
     */

    static final int INFO = 1;

    // #define CORRECTION pour faire une adaptation de la fo

    /*
     * on coupe si on a pas XX% de la dyn. de
     * l'energie
     */
    static final int SEUIL_NRJ = 40;
    /*
     * on coupe si on a plus que YY% de la dyn
     * de dpz
     */
    static final int SEUIL_DPZ = 50;

    /* seuil pour le calcul de la dpz */
    static final int EPSILON = 10;

    static final int SEUIL_VOISE = 7;

    static final int AVANCE = 2;

    /* max. number of amdf values per trame */
    static final int cst_pics_amdf = 5;

    /*
     * percentage of pitch-variability
     * between 2 values
     */
    static final int PITCH_VARIABILITY = 25;

    /*
     * number of time it should apply the
     * high frequency filter
     */
    static final int FILTRE_PASSE_BAS = 5;

    static final int POURCENT(int n, int t) {
        return t != 0 ? n * 100 / t : 0;
    }

    static final double CARRE(double a) {
        return a * a;
    }

    static final int TO_FREQ(int p) {
        return p != 0 ? cst_freq_ech / p : 0;
    }

    static final int SEEK_SET = 0;
    static final int SEEK_END = 2;

    static final int MAX_ENTIER = 2147483;

    final boolean NON_VOISEE(int t) {
        return Vois[t] < SEUIL_VOISE;
    }

    final boolean VOISEE(int t) {
        return Vois[t] >= SEUIL_VOISE;
    }

    /*
     * Quelques declarations anodines
     */

    class RESULT {
        int total;
        int rang;
    }

    class ZONE {
        int debut, fin, ancrage;
        ZONE suiv, pred;
    }

    /*
     * LES VARIABLES GLOBALES DU DETECTEUR
     */

    private int min_fo, max_fo, min_nrj, max_nrj, min_dpz, max_dpz;
    private int nb_coupe = 0, debug, quick;
    private int seuil_dpz, seuil_nrj;

    /*
     * 
     * hamming window length = 2.5 frequency / min_fo (16kHz,60Hz = 667) window
     * skip = frequency / 100 (16kHz = 160)
     */

    static int cst_freq_ech, cst_freq_coupure, cst_length_hamming, cst_step_hamming, cst_point_par_trame, cst_step_min, cst_step_max;

    private RESULT[][] Coeff_Amdf = new RESULT[cst_pics_amdf][];

    private ZONE zone;

    private double[] Hamming;

    private int max_amdf, min_amdf, amplitude_amdf;

    private short[] Nrj, Dpz, Vois, Fo;

    private float[] Signal;

    private int[][] Resultat;

    /** */
    private void init(int frequence, int cst_min_fo, int cst_max_fo) {
        cst_freq_coupure = 800;
        cst_freq_ech = frequence;
        cst_length_hamming = ((int) (2.5 * cst_freq_ech) / cst_min_fo);
        cst_step_hamming = cst_point_par_trame = (cst_freq_ech / 100);
        cst_step_min = (cst_freq_ech / cst_max_fo);
        cst_step_max = (cst_freq_ech / cst_min_fo);

        if (debug > 1)
            System.err.printf("sampling:%d, hamming size=%d, hamming overlap=%d\n", cst_freq_ech, cst_length_hamming, cst_step_hamming);
    }

    /** */
//    private void filtre_passe_bas(int frequence, int longueur) {
//        int i;
//        double coeff, delai = 0.0;
//        coeff = (Math.PI * 2 * frequence) / cst_freq_ech;
//        for (i = 0; i < longueur; i++) {
//            Signal[i] = (short) (delai = (Signal[i] * coeff) + (delai * (1 - coeff)));
//        }
//    }

    /** */
    private void libere_zone(ZONE zone) {
        ZONE l;

        while (zone != null) {
            l = zone.suiv;
            zone = l;
        }
    }

    /** */
    private void libere_coeff_amdf() {
    }

    /** */
    private int voisement_par_profondeur_des_pics(int imin, int[] result, int length) {
        int gauche, droite, i;

        // emergence gauche
        for (i = imin; i > 0 && result[i] <= result[i - 1]; i--) {
        }
        gauche = result[i] - result[imin];
        // emergence droite
        for (i = imin; i < length - 1 && result[i] <= result[i + 1]; i++) {
        }
        droite = result[i] - result[imin];

        return Math.min(droite, gauche);
    }

    /** */
    private void ranger_minimum(int trame, int valeur, int rang) {
        int i, j;

        for (i = 0; (i < cst_pics_amdf) && (valeur >= Coeff_Amdf[i][trame].total); i++) {
        }
        if (i < cst_pics_amdf) {
            for (j = cst_pics_amdf - 1; j > i; j--) {
                Coeff_Amdf[j][trame] = Coeff_Amdf[j - 1][trame];
            }
            Coeff_Amdf[i][trame].total = valeur;
            Coeff_Amdf[i][trame].rang = rang;
        }
    }

    /** */
    static final int homothetie(int y, int minval, int amp) {
        return amp != 0 ? y - minval * 200 / amp : 0;
    }

    private void retiens_n_pics(int[] result, int trame, int length, int[] profondeur) {

        int i, prof, profondeur_locale = 0;
        int[] result_local;
        int minVal, maxVal;

        // init
        for (i = 0; i < cst_pics_amdf; i++) {
            Coeff_Amdf[i][trame].total = MAX_ENTIER;
            Coeff_Amdf[i][trame].rang = -1;
        }

        // minimum
        for (i = 0; i < length - 1;) {
            while ((i < length - 1) && (result[i] > result[i + 1]))
                i++; /* on descend */
            if (i != 0 && i < length - 1)
                ranger_minimum(trame, result[i], i + cst_step_min);
            while ((i < length - 1) && (result[i] <= result[i + 1]))
                i++; /* on monte */
        }

        // Recherche du voisement par profondeur du pic minimum

        prof = profondeur[0] = 0;

        {
            result_local = new int[length];
            maxVal = 0;
            minVal = MAX_ENTIER;
            for (i = 0; i < length; i++) {
                if (result[i] > maxVal) {
                    maxVal = result[i];
                }
                if (result[i] < minVal) {
                    minVal = result[i];
                }
            }
            if (debug > 1) {
                System.err.printf("DYN AMDF[%d] : %d - %d (%d)  ", trame, minVal, maxVal, maxVal - minVal);
            }
        }

        for (i = 0; i < length; i++) {
            result_local[i] = homothetie(result[i], minVal, maxVal - minVal);
            result[i] = homothetie(result[i], min_amdf, amplitude_amdf);
        }

        for (i = 0; i < cst_pics_amdf; i++)
            if (Coeff_Amdf[i][trame].rang != -1) {
                prof = voisement_par_profondeur_des_pics(Coeff_Amdf[i][trame].rang - cst_step_min, result, length);
                if (prof > profondeur[0]) {
                    profondeur[0] = prof;
                }
                prof = voisement_par_profondeur_des_pics(Coeff_Amdf[i][trame].rang - cst_step_min, result_local, length);
                if (prof > profondeur_locale) {
                    profondeur_locale = prof;
                }
            }

        Vois[trame] = (short) profondeur_locale;
        if (debug > 1) {
            System.err.printf("---. %d\n", profondeur[0]);
        }
// #undef homothetie
    }

    /* *********************************************************************** */

    private void calcul_voisement(int nb_trames) {
        int trame, length;
        int[] profondeur = new int[1];

        amplitude_amdf = max_amdf - min_amdf;

        length = cst_step_max - cst_step_min + 1;

        for (trame = 0; trame < nb_trames; trame++) {
            if (quick != 0 && (Nrj[trame] < seuil_nrj) && (Dpz[trame] > seuil_dpz)) {
                Vois[trame] = 0;
            } else {
                retiens_n_pics(Resultat[trame], trame, length, profondeur);
                Vois[trame] = (short) profondeur[0];
            }
        }
    }

    /** */
    private void precalcul_hamming() {
        double pas;
        int i;

        pas = Math.PI * 2 / cst_length_hamming;
        for (i = 0; i < cst_length_hamming; i++) {
            Hamming[i] = 0.54 - 0.46 * Math.cos(i * pas);
        }
    }

    /** */
    private void amdf(Sound s, int i, int[] Hammer, int[] result, int nrj, int start) {
        final int FACTEUR = 50;

        int decal, j, l, somme, k;
        double coeff, delai;
        /* static */double[] odelai = new double[FILTRE_PASSE_BAS];

        Snack_GetSoundData(s, start + i, Signal, cst_length_hamming);

        if (i == 0) {
            for (k = 0; k < FILTRE_PASSE_BAS; k++) {
                odelai[k] = 0.0;
            }
        }

        for (k = 0; k < FILTRE_PASSE_BAS; k++) {
            delai = odelai[k];
            coeff = (Math.PI * 2 * cst_freq_coupure) / cst_freq_ech;
            for (j = 0; j < cst_length_hamming; j++) {
                Signal[j] = (float) (delai = Signal[j] * coeff + (delai * (1 - coeff)));
            }
            odelai[k] = Signal[cst_step_hamming - 1];
        }

        for (j = 0; j < cst_length_hamming; j++) {
            Hammer[j] = (int) (Hamming[j] * Signal[j]);
        }

        // la double boucle couteuse
        for (decal = cst_step_min; decal <= cst_step_max; decal++) {
            for (somme = l = 0, j = decal; j < cst_length_hamming; j++, l++) {
                somme += Math.abs(Hammer[j] - Hammer[l]);
            }
            result[decal - cst_step_min] = ((somme * FACTEUR) / ((cst_length_hamming - decal)));
        }
// #undef FACTEUR
    }

    /** */
    private int parametre_amdf(Sound s, Tcl.Interp interp, int start, int longueur, int[] nb_trames, int[] Hammer) {
        int j, i, length, trame;

        max_amdf = 0;
        min_amdf = MAX_ENTIER;

        length = cst_step_max - cst_step_min + 1;

        for (trame = i = 0; i < longueur; i += cst_step_hamming, trame++) {
            if (i > s.length - cst_length_hamming)
                break;
            if (i > longueur - cst_length_hamming / 2)
                break;
            if (quick != 0 && (Nrj[trame] < seuil_nrj) && (Dpz[trame] > seuil_dpz)) {
            } else {
                amdf(s, i, Hammer, Resultat[trame], (Nrj[trame] != 0) ? Nrj[trame] : 1, start);
                for (j = 0; j < length; j++) {
                    if (Resultat[trame][j] > max_amdf) {
                        max_amdf = Resultat[trame][j];
                    }
                    if (Resultat[trame][j] < min_amdf) {
                        min_amdf = Resultat[trame][j];
                    }
                }
            }
            if ((trame % 20) == 19) {
                int res = Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Computing pitch", 0.05 + 0.95 * i / longueur);
                if (res != 0) {
                    return -1;
                }
            }
        }
        Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Computing pitch", 1.0);
        nb_trames[0] = trame;

        if (debug != 0) {
            System.err.printf("min_amdf=%d, max_amdf=%d\n", min_amdf, max_amdf);
        }
        return 0;
    }

    /**
     * LA MAGOUILLE DE LA RECHERCHE DU MEILLEUR CHEMIN
     */
    private int interpolation(int x1, int y1, int x2, int y2, int x) {
        int a, b, delta;

        /* y = ax+b */

        if (x == x1) {
            return y1;
        }
        if (x == x2) {
            return y2;
        }

        if ((delta = x1 - x2) != 0) {
            a = y1 - y2;
            b = x1 * y2 - x2 * y1;

            return (a * x + b) / delta;
        }
        return 0;
    }

    /** */
    private boolean ecart_correct(int avant, int apres) {
        int a, b;

        a = cst_freq_ech / Math.max(avant, apres);
        b = cst_freq_ech / Math.min(avant, apres);

        return b <= (a * (100 + PITCH_VARIABILITY)) / 100;
    }

    private static final int TABLEAU2(int t, int reference, RESULT[] table) {
        return Math.abs(table[t].rang - reference);
    }

    /**
     * selection du pic le plus proche d'une reference
     */
    private void trier(int trame, int reference, RESULT[] table) {
        int t, bulle = 1;

        for (t = 0; t < cst_pics_amdf; t++) {
            table[t] = Coeff_Amdf[t][trame];
        }

        while (bulle != 0) {
            for (bulle = t = 0; t < cst_pics_amdf - 1; t++) {
                if ((table[t].rang == -1 && table[t + 1].rang != -1) || (TABLEAU2(t, reference, table) > TABLEAU2(t + 1, reference, table) && table[t + 1].rang != -1)) {
                    RESULT temp;
                    temp = table[t + 1];
                    table[t + 1] = table[t];
                    table[t] = temp;
                    bulle = 1;
                }
            }
        }
// #undef TABLEAU
    }

    /* */

    private static final int SEUIL = 20;

//    private static final int MAUVAIS_CHEMIN() {
//        return 0;
//    }

    private static final boolean ECART_CORRECT(int a, int b) {
        final int POURCENTAGE = 25;
        return (Math.max(a, b) - Math.min(a, b)) < Math.min(SEUIL, (POURCENTAGE * Math.min(a, b) / 100));
    }

    private void extension_zone(int accrochage, int debut, int fin, int to, short[] tableau) {

        int trame, avant, i, m, j;
        RESULT[][] table = new RESULT[AVANCE + 1][cst_pics_amdf];
        RESULT[] normal = new RESULT[cst_pics_amdf];

        tableau[accrochage] = (short) (avant = to);

        // on etend a gauche

        for (trame = accrochage; trame > debut;) {
            table[0][0].rang = avant;
            for (i = 1; i <= AVANCE && trame - i >= debut; i++) {
                trier(trame - i, table[i - 1][0].rang, table[i]);
            }
            trier(m = Math.max(trame - AVANCE, debut), avant, normal);

            i--;

            if (table[i][0].rang != -1 && normal[0].rang != -1 && table[i][0].rang != 0 && normal[0].rang != 0 && avant != 0 && Math.abs(TO_FREQ(table[i][0].rang) - TO_FREQ(avant)) > Math.abs(TO_FREQ(normal[0].rang) - TO_FREQ(avant))) {
                for (j = 1; j <= i; j++) {
                    tableau[trame - j] = (short) interpolation(trame, avant, m, normal[0].rang, trame - j);
                }
                trame -= i;
                avant = normal[0].rang;
            } else {
                if (table[1][0].rang != 0 && ECART_CORRECT(TO_FREQ(avant), TO_FREQ(table[1][0].rang))) {
                    avant = tableau[--trame] = (short) table[1][0].rang;
                } else {
                    tableau[--trame] = 0;
                }
            }
        }

        // on etend a droite

        avant = tableau[accrochage];
        for (trame = accrochage; trame < fin;) {
            table[0][0].rang = avant;
            for (i = 1; i <= AVANCE && trame + i <= fin; i++) {
                trier(trame + i, table[i - 1][0].rang, table[i]);
            }
            trier(m = Math.min(trame + AVANCE, fin), avant, normal);

            i--;

            if (table[i][0].rang != -1 && normal[0].rang != -1 && table[i][0].rang != 0 && normal[0].rang != 0 && avant != 0 && Math.abs(TO_FREQ(table[i][0].rang) - TO_FREQ(avant)) > Math.abs(TO_FREQ(normal[0].rang) - TO_FREQ(avant))) {
                for (j = 1; j <= i; j++) {
                    tableau[trame + j] = (short) interpolation(trame, avant, m, normal[0].rang, trame + j);
                }
                trame += i;

                avant = normal[0].rang;
            } else {
                if (table[1][0].rang != 0 && ECART_CORRECT(TO_FREQ(avant), TO_FREQ(table[1][0].rang))) {
                    avant = tableau[++trame] = (short) table[1][0].rang;
                } else {
                    tableau[++trame] = 0;
                }
            }
        }
// #undef MAUVAIS_CHEMIN
// #undef SEUIL
// #undef ECART_CORRECT
// #undef POURCENTAGE
    }

    /** */
    private void recupere_trou(ZONE l, short[] tableau) {
        int trame, avant;

        for (trame = l.debut; trame < l.fin;) {
            for (; trame <= l.fin && tableau[trame] != 0; trame++) {
            }
            for (avant = trame; trame <= l.fin && tableau[trame] == 0; trame++) {
            }
            while (avant > l.debut && trame < l.fin && (tableau[trame] == 0 || !ecart_correct(tableau[avant - 1], tableau[trame]))) {
                trame++;
            }
            if (avant > l.debut && trame < l.fin && ecart_correct(tableau[avant - 1], tableau[trame])) {
                int debut, fin, i;

                debut = avant - 1;
                fin = trame;
                for (i = debut + 1; i < fin; i++) {
                    tableau[i] = (short) interpolation(debut, tableau[debut], fin, tableau[fin], i);
                }
            }
        }

    }

    /** */
    private static final boolean PROCHE_MOYENNE(int e, int bande, int To_Moyenne) {
        return e >= To_Moyenne - bande && e <= To_Moyenne + bande;
    }

    private int point_accrochage(int debut, int fin, int[] accrochage, int To_Moyenne) {
        final int POURCENTAGE = 30;

        int bande, valeur, trame, pic;

        bande = To_Moyenne * POURCENTAGE / 100;
        valeur = MAX_ENTIER;

        for (pic = 0; pic <= 1; pic++) {
            for (trame = debut; trame <= fin; trame++) {
                if (Math.abs(Coeff_Amdf[pic][trame].rang - To_Moyenne) < Math.abs(valeur - To_Moyenne)) {
                    valeur = Coeff_Amdf[pic][trame].rang;
                    accrochage[0] = trame;
                }
            }
        }

        if (PROCHE_MOYENNE(valeur, bande, To_Moyenne)) {
            return valeur;
        }

        return 0;
// #undef POURCENTAGE
// #undef PROCHE_MOYENNE
    }

    /** */
    private int cherche_chemin(ZONE l, int To_moyen) {
        int[] accrochage = new int[1];

        if ((l.ancrage = point_accrochage(l.debut, l.fin, accrochage, To_moyen)) != 0) {
            return 0;
        }

        extension_zone(accrochage[0], l.debut, l.fin, l.ancrage, Fo);
        recupere_trou(l, Fo);
        return 1;
    }

    /** */
    private void calcul_courbe_fo(int nb_trames, int[] To_Moyen) {
        int t, memo, debut;
        ZONE l;

        // on met tout a zero au debut

        memo = To_Moyen[0];

        for (debut = 0; debut < nb_trames; debut++) {
            Fo[debut] = 0;
        }

        for (l = zone; l != null; l = l.suiv) {
            if (cherche_chemin(l, To_Moyen[0]) == 0) {
                for (t = l.debut; t <= l.fin; t++) {
                    Fo[t] = 0;
                }
            } else {
// #ifdef CORRECTION
                int cum = 0, nb = 0;

                for (t = l.fin; (t >= l.debut); t--) {
                    cum += Fo[t];
                    nb++;
                }
                To_Moyen[0] = ((2 * cum) + (nb * memo)) / (3 * nb);
                if (debug != 0) {
                    System.err.printf("correction moyenne : %d\n", To_Moyen[0]);
                }
// #endif
            }
        }

        // on recalcule la moyenne sur le chemin choisi

        min_fo = MAX_ENTIER;
        max_fo = 0;

        for ((To_Moyen[0]) = debut = t = 0; t < nb_trames; t++) {
            if (Fo[t] != 0) {
                (To_Moyen[0]) += Fo[t];
                Fo[t] = (short) (cst_freq_ech / Fo[t]);
                if (max_fo < Fo[t]) {
                    max_fo = Fo[t];
                }
                if (min_fo > Fo[t]) {
                    min_fo = Fo[t];
                }
                debut++;
            }
        }
        if (debut != 0) {
            (To_Moyen[0]) /= debut;
        }

        if (debug != 0) {
            System.err.printf("MOYENNE FINALE : %d (fo=%d)\n", To_Moyen[0], cst_freq_ech / To_Moyen[0]);
        }
    }

    /** */
    private static final int POURCENTAGE = 30;

//    private static final boolean ACCEPTABLE(int valeur, int moyenne, int pourcentage) {
//        return valeur > moyenne - pourcentage && valeur < moyenne + pourcentage;
//    }

    private static final int TABLEAU3(int t, RESULT[] table, int To_Moyenne) {
        return Math.abs(table[t].rang - (To_Moyenne));
    }

    private void calcul_fo_moyen(int nb_trames, int[] To_Moyenne) {
        int trame, nb, bulle;
        RESULT[] table;
        int To_Moyenne_Corrigee;

        table = new RESULT[nb_trames];

        for (To_Moyenne[0] = trame = nb = 0; trame < nb_trames; trame++) {
            if (VOISEE(trame)) {
                table[nb++] = Coeff_Amdf[0][trame];
                To_Moyenne[0] += Coeff_Amdf[0][trame].rang;
            }
        }

        To_Moyenne[0] = (nb != 0) ? ((To_Moyenne[0]) / nb) : 1;

        if (debug != 0) {
            System.err.printf("To moyen non corrige : %d (fo=%d) \n", To_Moyenne[0], cst_freq_ech / To_Moyenne[0]);
        }

        // correction de la valeur de fo

        for (bulle = 1; bulle != 0;) {
            for (bulle = trame = 0; trame < nb - 1; trame++) {
                if (TABLEAU3(trame, table, To_Moyenne[0]) > TABLEAU3(trame + 1, table, To_Moyenne[0])) {
                    RESULT temp;
                    bulle = 1;
                    temp = table[trame];
                    table[trame] = table[trame + 1];
                    table[trame + 1] = temp;
                }
            }
        }

        nb -= (POURCENTAGE * nb) / 100;

        for (To_Moyenne_Corrigee = 0, trame = 0; trame < nb; trame++) {
            To_Moyenne_Corrigee += table[trame].rang;
        }

        To_Moyenne_Corrigee = (nb != 0) ? (To_Moyenne_Corrigee / nb) : 1;

        // resultats

        To_Moyenne[0] = To_Moyenne_Corrigee;

        if (debug != 0) {
            System.err.printf("moyenne (a %d%% presque partout): %d (fo=%d)\n", 100 - POURCENTAGE, To_Moyenne[0], cst_freq_ech / To_Moyenne[0]);
        }
// #undef POURCENTAGE
// #undef ACCEPTABLE
// #undef TABLEAU
    }

    /** */
    private ZONE calcul_zones_voisees(int nb_trames) {
        int trame, debut;
        ZONE z = null;

        for (trame = 0; trame < nb_trames;) {
            while ((trame < nb_trames) && (NON_VOISEE(trame))) {
                trame++;
            }
            for (debut = trame; trame < nb_trames && VOISEE(trame); trame++) {
            }

            if (trame <= nb_trames && debut < trame) {
                ZONE t, l, pred;

                t = new ZONE();
                t.debut = debut;
                t.fin = trame - 1;
                t.ancrage = 0;
                t.suiv = null;

                for (l = z, pred = null; l != null; pred = l, l = l.suiv) {
                }

                t.pred = pred;
                if (pred != null) {
                    pred.suiv = t;
                } else {
                    z = t;
                }
            }
        }

        return z;
    }

    /** */
    private int calcul_nrj_dpz(Sound s, Tcl.Interp interp, int start, int longueur) {
        int J, JJ, m, i, j, trame, dpz;
        boolean sens;
        double nrj;

        max_nrj = max_dpz = 0;
        min_nrj = min_dpz = MAX_ENTIER;

        Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Computing pitch", 0.0);

        for (trame = i = 0; i < longueur; i += cst_step_hamming, trame++) {
            J = Math.min(s.length, (i + cst_length_hamming));
            JJ = J - 1;
            if (s.length < i + start + cst_length_hamming) {
                Snack_GetSoundData(s, i + start, Signal, s.length - i + start);
                for (j = s.length - i + start; j < cst_length_hamming; j++) {
                    Signal[j] = 0.0f;
                }
            } else {
                Snack_GetSoundData(s, i + start, Signal, cst_length_hamming);
            }

            /* ---- nrj ---- */
            for (nrj = 0.0, j = 0; j < J - i; j++)
                nrj += CARRE(Signal[j]);

            m = Nrj[trame] = (short) (10 * Math.log10(nrj));

            if (m > max_nrj) {
                max_nrj = m;
            }
            if (m < min_nrj) {
                min_nrj = m;
            }

            /* ---- dpz ---- */
            for (dpz = 0, j = 0; j < J - i; j++) {
                while ((j < J - i) && (Math.abs((int) Signal[j]) > EPSILON)) {
                    j++; /* looking just close to zero values */
                }
                if (j < J - i) {
                    dpz++;
                }
                sens = (((j - 1) >= 0) && (Signal[j - 1] > Signal[j]));
                if (sens) {
                    while ((j < JJ - i) && (Signal[j] > Signal[j + 1])) {
                        j++;
                    }
                } else {
                    while ((j < JJ - i) && (Signal[j] <= Signal[j + 1])) {
                        j++;
                    }
                }
            }
            m = Dpz[trame] = (short) dpz;

            if (m > max_dpz) {
                max_dpz = m;
            }
            if (m < min_dpz) {
                min_dpz = m;
            }

            if ((trame % 300) == 299) {
                int res = Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Computing pitch", 0.05 * i / longueur);
                if (res != 0) {
                    return -1;
                }
            }
        }

        seuil_nrj = min_nrj + (SEUIL_NRJ * (max_nrj - min_nrj)) / 100;
        seuil_dpz = min_dpz + (SEUIL_DPZ * (max_dpz - min_dpz)) / 100;

        if (debug != 0) {
            System.err.printf("dpz <%d,%d>, nrj <%d,%d> => Seuil nrj: %d, Seuil dpz: %d\n", min_dpz, max_dpz, min_nrj, max_nrj, seuil_nrj, seuil_dpz);
        }

        return trame;
    }

    /** */
//    private void adjust_signal(int longueur, int freq, int debug) {
//        int i;
//        long moy = 0, m;
//
//        for (i = 0; i < (m = Math.min(freq, longueur)); moy += Signal[i++]) {
//        }
//        moy /= (m != 0) ? m : 1;
//        if (debug != 0) {
//            System.err.printf("ajustement de %d dans tout le signal\n", moy);
//        }
//        if (moy != 0) {
//            for (i = 0; i < longueur; Signal[i++] -= (short) moy) {
//            }
//        }
//    }

    /** */
    private enum subOptions {
        START,
        END,
        F0MAX,
        F0MIN,
        PROGRESS,
        FRAME,
        METHOD,
        WINLEN
    }

    GetF0 getf0 = new GetF0();

    public int exec(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        int longueur, nb_trames, nb_trames_alloc;
        int[] To_Moyen = new int[1];
        int[] Hammer;
        int i;
        int fmin = 60, fmax = 400, lquick = 1/* , adjust = 0 */;
        int[] nbframes = new int[1];
        short minnrj, maxnrj, minfo, maxfo, mindpz, maxdpz;
        int arg, startpos = 0, endpos = -1, result, start;
        Tcl.Obj list;
        String[] subOptionStrings = {
            "-start", "-end", "-maxpitch", "-minpitch", "-progress", "-framelength", "-method", "-windowlength", null
        };

        if (s.debug > 0) {
            System.err.printf("Enter pitchCmd\n");
        }

        if (s.nchannels[0] != 1) {
            Tcl.AppendResult(interp, "pitch only works with Mono sounds", (String) null);
            return -1;
        }

        for (arg = 2; arg < objc; arg += 2) {
            String opt = Tcl.GetStringFromObj(objv[arg], null);
            String val = Tcl.GetStringFromObj(objv[arg + 1], null);

            if ("-method".equals(opt) && "esps".equalsIgnoreCase(val)) {
                getf0.Get_f0(s, interp, objc, objv);
                return 0;
            }
        }

        if (s.cmdPtr != null) {
            Tcl.DecrRefCount(s.cmdPtr);
            s.cmdPtr = null;
        }

        for (arg = 2; arg < objc; arg += 2) {
            int index;

            if (Tcl.GetIndexFromObj(interp, objv[arg], subOptionStrings, "option", 0, index) != 0) {
                return -1;
            }

            if (arg + 1 == objc) {
                Tcl.AppendResult(interp, "No argument given for ", subOptionStrings[index], " option", (String) null);
                return -1;
            }

            switch (subOptions.values()[index]) {
            case START: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], startpos) != 0)
                    return -1;
                break;
            }
            case END: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], endpos) != 0)
                    return -1;
                break;
            }
            case F0MAX: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], fmax) != 0)
                    return -1;
                if (fmax <= 50) {
                    Tcl.AppendResult(interp, "Maximum pitch must be > 50", null);
                    return -1;
                }
                break;
            }
            case F0MIN: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], fmin) != 0)
                    return -1;
                if (fmin <= 50) {
                    Tcl.AppendResult(interp, "Minimum pitch must be > 50", null);
                    return -1;
                }
                break;
            }
            case PROGRESS: {
                String str = Tcl.GetStringFromObj(objv[arg + 1], null);

                if (str.length() > 0) {
                    Tcl.IncrRefCount(objv[arg + 1]);
                    s.cmdPtr = objv[arg + 1];
                }
                break;
            }
            case FRAME: {
                break;
            }
            case METHOD: {
                break;
            }
            case WINLEN: {
                break;
            }
            }
        }
        if (fmax <= fmin) {
            Tcl.AppendResult(interp, "maxpitch must be > minpitch", null);
            return -1;
        }
        if (startpos < 0) {
            startpos = 0;
        }
        if (endpos >= (s.length - 1) || endpos == -1) {
            endpos = s.length - 1;
        }
        if (startpos > endpos) {
            return 0;
        }

        quick = lquick;
        init(s.samprate, fmin, fmax);

        start = startpos - (cst_length_hamming / 2);
        if (start < 0) {
            start = 0;
        }
        if (endpos - start + 1 < cst_length_hamming) {
            endpos = cst_length_hamming + start - 1;
            if (endpos >= s.length) {
                return 0;
            }
        }
        longueur = endpos - start + 1;

        Signal = new float[cst_length_hamming];

//        if (adjust) {
//            adjust_signal(longueur, s.samprate, s.debug);
//        }

        nb_trames = (longueur / cst_step_hamming) + 10;
        Nrj = new short[nb_trames];
        Dpz = new short[nb_trames];
        Vois = new short[nb_trames];
        Fo = new short[nb_trames];
        Resultat = new int[nb_trames][];

        for (i = 0; i < nb_trames; i++) {
            Resultat[i] = new int[cst_step_max - cst_step_min + 1];
        }
        nb_trames_alloc = nb_trames;
        nb_trames = nbframes[0] = calcul_nrj_dpz(s, interp, start, longueur);

        Hamming = new double[cst_length_hamming];
        Hammer = new int[cst_length_hamming];
        for (i = 0; i < cst_pics_amdf; i++) {
            Coeff_Amdf[i] = new RESULT[nb_trames];
        }

        precalcul_hamming();

        result = parametre_amdf(s, interp, start, longueur, nbframes, Hammer);

        if (result == 0) {
            if (debug != 0) {
                System.err.printf("nbframes=%d\n", nbframes);
            }
            calcul_voisement(nbframes[0]);
            zone = calcul_zones_voisees(nbframes[0]);
            calcul_fo_moyen(nbframes[0], To_Moyen);
            calcul_courbe_fo(nbframes[0], To_Moyen);

            minfo = (short) min_fo;
            maxfo = (short) max_fo;
            minnrj = (short) min_nrj;
            maxnrj = (short) max_nrj;
            mindpz = (short) min_dpz;
            maxdpz = (short) max_dpz;

            if (debug != 0 && quick != 0) {
                System.err.printf("%d trames coupees sur %d . %d %% (seuil nrj = %d, seuil dpz = %d) \n", nb_coupe, nbframes, POURCENT(nb_coupe, nbframes[0]), seuil_nrj, seuil_dpz);
            }
            libere_zone(zone);
        }
        libere_coeff_amdf();
        if (result == 1) {
            int n = cst_length_hamming / (2 * cst_step_hamming) - startpos / cst_step_hamming;

            list = Tcl.NewListObj(0, null);
            for (i = 0; i < n; i++) {
                Tcl.ListObjAppendElement(interp, list, Tcl.NewDoubleObj(0.0));
            }
            for (i = 0; i < nbframes[0]; i++) {
                Tcl.ListObjAppendElement(interp, list, Tcl.NewDoubleObj(Fo[i]));
            }
            Tcl.SetObjResult(interp, list);
        }

        if (s.debug > 0) {
            System.err.printf("Exit pitchCmd\n");
        }

        return 0;
    }

    private int cPitch(Sound s, Tcl.Interp interp, int[][] outlist, int[] length) {
        int longueur, nb_trames;
        int[] To_Moyen = new int[1];
        int[] Hammer;
        int i;
        int fmin = 60, fmax = 400, lquick = 1/* , adjust = 0 */;
        int[] nbframes = new int[1];
        short minnrj, maxnrj, minfo, maxfo, mindpz, maxdpz;
        int startpos = 0, endpos = -1, result, start;

        if (s.debug > 0) {
            System.err.printf("Enter pitchCmd\n");
        }

        if (fmax <= fmin) {
            Tcl.AppendResult(interp, "maxpitch must be > minpitch", null);
            return -1;
        }
        if (startpos < 0) {
            startpos = 0;
        }
        if (endpos >= (s.length - 1) || endpos == -1) {
            endpos = s.length - 1;
        }
        if (startpos > endpos) {
            return 0;
        }

        quick = lquick;
        init(s.samprate, fmin, fmax);

        start = startpos - (cst_length_hamming / 2);
        if (start < 0) {
            start = 0;
        }
        longueur = endpos - start + 1;

        Signal = new float[cst_length_hamming];

        nb_trames = (longueur / cst_step_hamming) + 10;
        Nrj = new short[nb_trames];
        Dpz = new short[nb_trames];
        Vois = new short[nb_trames];
        Fo = new short[nb_trames];
        Resultat = new int[nb_trames][];

        for (i = 0; i < nb_trames; i++) {
            Resultat[i] = new int[cst_step_max - cst_step_min + 1];
        }

        nb_trames = nbframes[0] = calcul_nrj_dpz(s, interp, start, longueur);

        Hamming = new double[cst_length_hamming];
        Hammer = new int[cst_length_hamming];
        for (i = 0; i < cst_pics_amdf; i++) {
            Coeff_Amdf[i] = new RESULT[nb_trames];
        }

        precalcul_hamming();

        result = parametre_amdf(s, interp, start, longueur, nbframes, Hammer);

        if (result == 0) {
            calcul_voisement(nbframes[0]);
            zone = calcul_zones_voisees(nbframes[0]);
            calcul_fo_moyen(nbframes[0], To_Moyen);
            calcul_courbe_fo(nbframes[0], To_Moyen);

            minfo = (short) min_fo;
            maxfo = (short) max_fo;
            minnrj = (short) min_nrj;
            maxnrj = (short) max_nrj;
            mindpz = (short) min_dpz;
            maxdpz = (short) max_dpz;

            libere_zone(zone);
        }
        libere_coeff_amdf();
        if (result == 0) {
            int n = cst_length_hamming / (2 * cst_step_hamming) - startpos / cst_step_hamming;
            int[] tmp = new int[n + nb_trames];
            for (i = 0; i < n; i++) {
                tmp[i] = 0;
            }
            for (i = n; i < nbframes[0] + n; i++) {
                tmp[i] = Fo[i - n];
            }
            outlist[0] = tmp;
            length[0] = nbframes[0] + n;
        }

        if (s.debug > 0) {
            System.err.printf("Exit pitchCmd\n");
        }

        return 0;
    }
}

/* */
