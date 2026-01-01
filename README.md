# Commande prédictive (MPC) pour la gestion énergétique d’un micro-réseau résidentiel (PV–Batterie–Réseau)

Ce projet présente une stratégie de **commande prédictive (MPC)** appliquée à la gestion énergétique d’un **micro-réseau résidentiel** composé d’une **production photovoltaïque (PV)**, d’une **batterie** et d’un **échange bidirectionnel avec le réseau**.  
L’objectif est de **minimiser le coût de l’énergie échangée avec le réseau** sous tarification **Heures Pleines / Heures Creuses (HP/HC)**, tout en respectant les **contraintes réseau** (abonnement) et les **contraintes de stockage** (SoC/énergie batterie).

---

## Objectifs du projet
- Modéliser un micro-réseau PV–batterie–réseau et ses flux de puissance.
- Implémenter une stratégie **MPC** basée sur un problème d’optimisation à horizon fini.
- Étudier la **sensibilité aux erreurs de prédiction** (charge et PV).
- Comparer des solveurs d’optimisation :  
  - **SQP** (local) via `fmincon`  
  - **PSO** (global) via `particleswarm`
- Réaliser une **étude techno-économique** de l’intérêt de la batterie en fonction de sa capacité (0–100 kWh) et du coût d’investissement.

---

## Architecture du dépôt
  .
  ├── code_matlab/
  │   ├── funobj.m                 # Fonction objectif (coût réseau + pénalisation contraintes molles)
  │   ├── systeme.m                # Modèle batterie + saturation + mise à jour Eb
  │   ├── script_1_cas_base.m       # 1) Cas de base (prédictions parfaites)
  │   ├── script_2_changer_donnees.m# 2) Étude avec données réelles vs données prédites
  │   ├── script_3_initialisation.m # 3) Effet du point initial (SQP)
  │   ├── script_4_PSO_vs_SQP.m     # 4) Comparaison PSO vs SQP (temps réel / performance)
  │   └── script_5_etude_batterie.m # 5) Étude de la capacité batterie + analyse économique
  │
  └── documentation/
      ├── rapport.pdf              # Rapport complet (modèles, MPC, résultats, analyses)
      └── fiche_resume.pdf         # Fiche A4 synthèse (1 page)


---

## Modèle résumé (idée générale)

- **Bilan de puissance**
  \[
  P_b = P_r + P_p - P_c
  \]
- **Dynamique énergie batterie**
  \[
  E_b(k+1)=E_b(k)+P_b(k)\frac{T_e}{60}
  \]
- **Contraintes**
  - Puissance réseau : \(P_{r,\min} \le P_r \le P_{r,\max}\)
  - Énergie batterie : \(0.1C_b \le E_b \le C_b\) (seuil bas pour limiter la dégradation)

- **MPC**
  - horizon typique : 24 h
  - pas : 10 min
  - optimisation à chaque pas, application de la première commande (receding horizon)

---

## Prérequis
- **MATLAB** (avec :
  - *Optimization Toolbox* pour `fmincon`
  - *Global Optimization Toolbox* pour `particleswarm`)

---

## Comment exécuter le projet (MATLAB)

1. Ouvrir MATLAB puis se placer dans le dossier :

2. Lancer les scripts dans l’ordre recommandé :

- **1) Cas de base**
- `script_1_cas_base.m`

- **2) Données réelles vs données prédites**
- `script_2_changer_donnees.m`

- **3) Effet de l’initialisation (SQP)**
- `script_3_initialisation.m`

- **4) Comparaison PSO vs SQP**
- `script_4_PSO_vs_SQP.m`

- **5) Étude de la batterie**
- `script_5_etude_batterie.m`

> Les figures et résultats d’analyse sont décrits en détail dans `documentation/rapport.pdf`.

---

## Fichiers clés (code)
- **`systeme.m`**
- Simule la dynamique de la batterie
- Applique les saturations (limites SoC/énergie)
- Retourne la puissance batterie réelle et l’état mis à jour

- **`funobj.m`**
- Calcule le coût total sur l’horizon MPC
- Intègre la pénalisation des contraintes molles sur la puissance réseau réelle

---

## Résultats (aperçu)
Le projet met en évidence :
- la capacité du MPC à effectuer un **arbitrage temporel** (charger/décharger selon PV + tarif),
- une **robustesse globale** aux erreurs de prédiction dans le cas étudié,
- l’intérêt d’un solveur global (**PSO**) pour réduire le coût, au prix d’un temps de calcul plus élevé,
- une **capacité utile** au-delà de laquelle le gain devient marginal, et une rentabilité dépendante du prix d’achat de la batterie.

Les résultats détaillés (courbes, tableaux, discussions) sont dans :
- `documentation/rapport.pdf`
- `documentation/fiche_resume.pdf`

---

## Auteur
**Mohamed IGADARNE** — Master 2 EEA (UPVD)

---

## Licence
Projet académique (usage pédagogique).
