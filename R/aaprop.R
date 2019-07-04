#' Normalized amino acids properties
#' 
#' Normalized (0-1) 554 amino acid properties as retreived from AAIndex database 
#' (release 9.1) enriched with contactivity of amino acids.
#' 
#' @name aaprop
#' @docType data
#' @format A data frames with 20 columns and 600 rows.
#' @details 
#' Following properties are included (AAIndex key: description of the property)
#' \describe{ 
#'  \item{ANDN920101}{alpha-CH chemical shifts (Andersen et al., 1992)}
#'  \item{ARGP820101}{Hydrophobicity index (Argos et al., 1982)}
#'  \item{ARGP820102}{Signal sequence helical potential (Argos et al., 1982)}
#'  \item{ARGP820103}{Membrane-buried preference parameters (Argos et al., 1982)}
#'  \item{BEGF750101}{Conformational parameter of inner helix (Beghin-Dirkx, 1975)}
#'  \item{BEGF750102}{Conformational parameter of beta-structure (Beghin-Dirkx,
#'  1975)}
#'  \item{BEGF750103}{Conformational parameter of beta-turn (Beghin-Dirkx, 1975)}
#'  \item{BHAR880101}{Average flexibility indices (Bhaskaran-Ponnuswamy, 1988)}
#'  \item{BIGC670101}{Residue volume (Bigelow, 1967)}
#'  \item{BIOV880101}{Information value for accessibility; average fraction 35\%
#'  (Biou et al., 1988)}
#'  \item{BIOV880102}{Information value for accessibility; average fraction 23\%
#'  (Biou et al., 1988)}
#'  \item{BROC820101}{Retention coefficient in TFA (Browne et al., 1982)}
#'  \item{BROC820102}{Retention coefficient in HFBA (Browne et al., 1982)}
#'  \item{BULH740101}{Transfer free energy to surface (Bull-Breese, 1974)}
#'  \item{BULH740102}{Apparent partial specific volume (Bull-Breese, 1974)}
#'  \item{BUNA790101}{alpha-NH chemical shifts (Bundi-Wuthrich, 1979)}
#'  \item{BUNA790102}{alpha-CH chemical shifts (Bundi-Wuthrich, 1979)}
#'  \item{BUNA790103}{Spin-spin coupling constants 3JHalpha-NH (Bundi-Wuthrich,
#'  1979)}
#'  \item{BURA740101}{Normalized frequency of alpha-helix (Burgess et al., 1974)}
#'  \item{BURA740102}{Normalized frequency of extended structure (Burgess et al.,
#'  1974)}
#'  \item{CHAM810101}{Steric parameter (Charton, 1981)}
#'  \item{CHAM820101}{Polarizability parameter (Charton-Charton, 1982)}
#'  \item{CHAM820102}{Free energy of solution in water, kcal/mole (Charton-Charton,
#'  1982)}
#'  \item{CHAM830101}{The Chou-Fasman parameter of the coil conformation
#'  (Charton-Charton, 1983)}
#'  \item{CHAM830102}{A parameter defined from the residuals obtained from the best
#'  correlation of the Chou-Fasman parameter of beta-sheet (Charton-Charton, 1983)}
#'  \item{CHAM830103}{The number of atoms in the side chain labelled 1+1
#'  (Charton-Charton, 1983)}
#'  \item{CHAM830104}{The number of atoms in the side chain labelled 2+1
#'  (Charton-Charton, 1983)}
#'  \item{CHAM830105}{The number of atoms in the side chain labelled 3+1
#'  (Charton-Charton, 1983)}
#'  \item{CHAM830106}{The number of bonds in the longest chain (Charton-Charton,
#'  1983)}
#'  \item{CHAM830107}{A parameter of charge transfer capability (Charton-Charton,
#'  1983)}
#'  \item{CHAM830108}{A parameter of charge transfer donor capability
#'  (Charton-Charton, 1983)}
#'  \item{CHOC750101}{Average volume of buried residue (Chothia, 1975)}
#'  \item{CHOC760101}{Residue accessible surface area in tripeptide (Chothia,
#'  1976)}
#'  \item{CHOC760102}{Residue accessible surface area in folded protein (Chothia,
#'  1976)}
#'  \item{CHOC760103}{Proportion of residues 95\% buried (Chothia, 1976)}
#'  \item{CHOC760104}{Proportion of residues 100\% buried (Chothia, 1976)}

#'  \item{CHOP780101}{Normalized frequency of beta-turn (Chou-Fasman, 1978a)}
#'  \item{CHOP780201}{Normalized frequency of alpha-helix (Chou-Fasman, 1978b)}
#'  \item{CHOP780202}{Normalized frequency of beta-sheet (Chou-Fasman, 1978b)}
#'  \item{CHOP780203}{Normalized frequency of beta-turn (Chou-Fasman, 1978b)}
#'  \item{CHOP780204}{Normalized frequency of N-terminal helix (Chou-Fasman,
#'  1978b)}
#'  \item{CHOP780205}{Normalized frequency of C-terminal helix (Chou-Fasman,
#'  1978b)}
#'  \item{CHOP780206}{Normalized frequency of N-terminal non helical region
#'  (Chou-Fasman, 1978b)}
#'  \item{CHOP780207}{Normalized frequency of C-terminal non helical region
#'  (Chou-Fasman, 1978b)}
#'  \item{CHOP780208}{Normalized frequency of N-terminal beta-sheet (Chou-Fasman,
#'  1978b)}
#'  \item{CHOP780209}{Normalized frequency of C-terminal beta-sheet (Chou-Fasman,
#'  1978b)}
#'  \item{CHOP780210}{Normalized frequency of N-terminal non beta region
#'  (Chou-Fasman, 1978b)}
#'  \item{CHOP780211}{Normalized frequency of C-terminal non beta region
#'  (Chou-Fasman, 1978b)}
#'  \item{CHOP780212}{Frequency of the 1st residue in turn (Chou-Fasman, 1978b)}
#'  \item{CHOP780213}{Frequency of the 2nd residue in turn (Chou-Fasman, 1978b)}
#'  \item{CHOP780214}{Frequency of the 3rd residue in turn (Chou-Fasman, 1978b)}
#'  \item{CHOP780215}{Frequency of the 4th residue in turn (Chou-Fasman, 1978b)}
#'  \item{CHOP780216}{Normalized frequency of the 2nd and 3rd residues in turn
#'  (Chou-Fasman, 1978b)}
#'  \item{CIDH920101}{Normalized hydrophobicity scales for alpha-proteins (Cid et
#'  al., 1992)}
#'  \item{CIDH920102}{Normalized hydrophobicity scales for beta-proteins (Cid et
#'  al., 1992)}
#'  \item{CIDH920103}{Normalized hydrophobicity scales for alpha+beta-proteins (Cid
#'  et al., 1992)}
#'  \item{CIDH920104}{Normalized hydrophobicity scales for alpha/beta-proteins (Cid
#'  et al., 1992)}
#'  \item{CIDH920105}{Normalized average hydrophobicity scales (Cid et al., 1992)}
#'  \item{COHE430101}{Partial specific volume (Cohn-Edsall, 1943)}
#'  \item{CRAJ730101}{Normalized frequency of middle helix (Crawford et al., 1973)}
#'  \item{CRAJ730102}{Normalized frequency of beta-sheet (Crawford et al., 1973)}
#'  \item{CRAJ730103}{Normalized frequency of turn (Crawford et al., 1973)}
#'  \item{DAWD720101}{Size (Dawson, 1972)}
#'  \item{DAYM780101}{Amino acid composition (Dayhoff et al., 1978a)}
#'  \item{DAYM780201}{Relative mutability (Dayhoff et al., 1978b)}
#'  \item{DESM900101}{Membrane preference for cytochrome b: MPH89 (Degli Esposti et
#'  al., 1990)}
#'  \item{DESM900102}{Average membrane preference: AMP07 (Degli Esposti et al.,
#'  1990)}
#'  \item{EISD840101}{Consensus normalized hydrophobicity scale (Eisenberg, 1984)}
#'  \item{EISD860101}{Solvation free energy (Eisenberg-McLachlan, 1986)}
#'  \item{EISD860102}{Atom-based hydrophobic moment (Eisenberg-McLachlan, 1986)}
#'  \item{EISD860103}{Direction of hydrophobic moment (Eisenberg-McLachlan, 1986)}
#'  \item{FASG760101}{Molecular weight (Fasman, 1976)}
#'  \item{FASG760102}{Melting point (Fasman, 1976)}
#'  \item{FASG760103}{Optical rotation (Fasman, 1976)}
#'  \item{FASG760104}{pK-N (Fasman, 1976)}
#'  \item{FASG760105}{pK-C (Fasman, 1976)}
#'  \item{FAUJ830101}{Hydrophobic parameter pi (Fauchere-Pliska, 1983)}
#'  \item{FAUJ880101}{Graph shape index (Fauchere et al., 1988)}
#'  \item{FAUJ880102}{Smoothed upsilon steric parameter (Fauchere et al., 1988)}
#'  \item{FAUJ880103}{Normalized van der Waals volume (Fauchere et al., 1988)}
#'  \item{FAUJ880104}{STERIMOL length of the side chain (Fauchere et al., 1988)}
#'  \item{FAUJ880105}{STERIMOL minimum width of the side chain (Fauchere et al.,
#'  1988)}
#'  \item{FAUJ880106}{STERIMOL maximum width of the side chain (Fauchere et al.,
#'  1988)}
#'  \item{FAUJ880107}{N.m.r. chemical shift of alpha-carbon (Fauchere et al.,
#'  1988)}
#'  \item{FAUJ880108}{Localized electrical effect (Fauchere et al., 1988)}
#'  \item{FAUJ880109}{Number of hydrogen bond donors (Fauchere et al., 1988)}
#'  \item{FAUJ880110}{Number of full nonbonding orbitals (Fauchere et al., 1988)}
#'  \item{FAUJ880111}{Positive charge (Fauchere et al., 1988)}
#'  \item{FAUJ880112}{Negative charge (Fauchere et al., 1988)}
#'  \item{FAUJ880113}{pK-a(RCOOH) (Fauchere et al., 1988)}
#'  \item{FINA770101}{Helix-coil equilibrium constant (Finkelstein-Ptitsyn, 1977)}
#'  \item{FINA910101}{Helix initiation parameter at posision i-1 (Finkelstein et
#'  al., 1991)}
#'  \item{FINA910102}{Helix initiation parameter at posision i,i+1,i+2 (Finkelstein
#'  et al., 1991)}
#'  \item{FINA910103}{Helix termination parameter at posision j-2,j-1,j
#'  (Finkelstein et al., 1991)}
#'  \item{FINA910104}{Helix termination parameter at posision j+1 (Finkelstein et
#'  al., 1991)}
#'  \item{GARJ730101}{Partition coefficient (Garel et al., 1973)}
#'  \item{GEIM800101}{Alpha-helix indices (Geisow-Roberts, 1980)}
#'  \item{GEIM800102}{Alpha-helix indices for alpha-proteins (Geisow-Roberts,
#'  1980)}
#'  \item{GEIM800103}{Alpha-helix indices for beta-proteins (Geisow-Roberts, 1980)}
#'  \item{GEIM800104}{Alpha-helix indices for alpha/beta-proteins (Geisow-Roberts,
#'  1980)}
#'  \item{GEIM800105}{Beta-strand indices (Geisow-Roberts, 1980)}
#'  \item{GEIM800106}{Beta-strand indices for beta-proteins (Geisow-Roberts, 1980)}
#'  \item{GEIM800107}{Beta-strand indices for alpha/beta-proteins (Geisow-Roberts,
#'  1980)}
#'  \item{GEIM800108}{Aperiodic indices (Geisow-Roberts, 1980)}
#'  \item{GEIM800109}{Aperiodic indices for alpha-proteins (Geisow-Roberts, 1980)}
#'  \item{GEIM800110}{Aperiodic indices for beta-proteins (Geisow-Roberts, 1980)}
#'  \item{GEIM800111}{Aperiodic indices for alpha/beta-proteins (Geisow-Roberts,
#'  1980)}
#'  \item{GOLD730101}{Hydrophobicity factor (Goldsack-Chalifoux, 1973)}
#'  \item{GOLD730102}{Residue volume (Goldsack-Chalifoux, 1973)}
#'  \item{GRAR740101}{Composition (Grantham, 1974)}
#'  \item{GRAR740102}{Polarity (Grantham, 1974)}
#'  \item{GRAR740103}{Volume (Grantham, 1974)}
#'  \item{GUYH850101}{Partition energy (Guy, 1985)}
#'  \item{HOPA770101}{Hydration number (Hopfinger, 1971), Cited by Charton-Charton
#'  (1982)}
#'  \item{HOPT810101}{Hydrophilicity value (Hopp-Woods, 1981)}
#'  \item{HUTJ700101}{Heat capacity (Hutchens, 1970)}
#'  \item{HUTJ700102}{Absolute entropy (Hutchens, 1970)}
#'  \item{HUTJ700103}{Entropy of formation (Hutchens, 1970)}
#'  \item{ISOY800101}{Normalized relative frequency of alpha-helix (Isogai et al.,
#'  1980)}
#'  \item{ISOY800102}{Normalized relative frequency of extended structure (Isogai
#'  et al., 1980)}
#'  \item{ISOY800103}{Normalized relative frequency of bend (Isogai et al., 1980)}
#'  \item{ISOY800104}{Normalized relative frequency of bend R (Isogai et al.,
#'  1980)}
#'  \item{ISOY800105}{Normalized relative frequency of bend S (Isogai et al.,
#'  1980)}
#'  \item{ISOY800106}{Normalized relative frequency of helix end (Isogai et al.,
#'  1980)}
#'  \item{ISOY800107}{Normalized relative frequency of double bend (Isogai et al.,
#'  1980)}
#'  \item{ISOY800108}{Normalized relative frequency of coil (Isogai et al., 1980)}
#'  \item{JANJ780101}{Average accessible surface area (Janin et al., 1978)}
#'  \item{JANJ780102}{Percentage of buried residues (Janin et al., 1978)}
#'  \item{JANJ780103}{Percentage of exposed residues (Janin et al., 1978)}
#'  \item{JANJ790101}{Ratio of buried and accessible molar fractions (Janin, 1979)}
#'  \item{JANJ790102}{Transfer free energy (Janin, 1979)}
#'  \item{JOND750101}{Hydrophobicity (Jones, 1975)}
#'  \item{JOND750102}{pK (-COOH) (Jones, 1975)}
#'  \item{JOND920101}{Relative frequency of occurrence (Jones et al., 1992)}
#'  \item{JOND920102}{Relative mutability (Jones et al., 1992)}
#'  \item{JUKT750101}{Amino acid distribution (Jukes et al., 1975)}
#'  \item{JUNJ780101}{Sequence frequency (Jungck, 1978)}
#'  \item{KANM800101}{Average relative probability of helix (Kanehisa-Tsong, 1980)}
#'  \item{KANM800102}{Average relative probability of beta-sheet (Kanehisa-Tsong,
#'  1980)}
#'  \item{KANM800103}{Average relative probability of inner helix (Kanehisa-Tsong,
#'  1980)}
#'  \item{KANM800104}{Average relative probability of inner beta-sheet
#'  (Kanehisa-Tsong, 1980)}
#'  \item{KARP850101}{Flexibility parameter for no rigid neighbors (Karplus-Schulz,
#'  1985)}
#'  \item{KARP850102}{Flexibility parameter for one rigid neighbor (Karplus-Schulz,
#'  1985)}
#'  \item{KARP850103}{Flexibility parameter for two rigid neighbors
#'  (Karplus-Schulz, 1985)}
#'  \item{KHAG800101}{The Kerr-constant increments (Khanarian-Moore, 1980)}
#'  \item{KLEP840101}{Net charge (Klein et al., 1984)}
#'  \item{KRIW710101}{Side chain interaction parameter (Krigbaum-Rubin, 1971)}
#'  \item{KRIW790101}{Side chain interaction parameter (Krigbaum-Komoriya, 1979)}
#'  \item{KRIW790102}{Fraction of site occupied by water (Krigbaum-Komoriya, 1979)}
#'  \item{KRIW790103}{Side chain volume (Krigbaum-Komoriya, 1979)}
#'  \item{KYTJ820101}{Hydropathy index (Kyte-Doolittle, 1982)}
#'  \item{LAWE840101}{Transfer free energy, CHP/water (Lawson et al., 1984)}
#'  \item{LEVM760101}{Hydrophobic parameter (Levitt, 1976)}
#'  \item{LEVM760102}{Distance between C-alpha and centroid of side chain (Levitt,
#'  1976)}
#'  \item{LEVM760103}{Side chain angle theta(AAR) (Levitt, 1976)}
#'  \item{LEVM760104}{Side chain torsion angle phi(AAAR) (Levitt, 1976)}
#'  \item{LEVM760105}{Radius of gyration of side chain (Levitt, 1976)}
#'  \item{LEVM760106}{van der Waals parameter R0 (Levitt, 1976)}
#'  \item{LEVM760107}{van der Waals parameter epsilon (Levitt, 1976)}
#'  \item{LEVM780101}{Normalized frequency of alpha-helix, with weights (Levitt,
#'  1978)}
#'  \item{LEVM780102}{Normalized frequency of beta-sheet, with weights (Levitt,
#'  1978)}
#'  \item{LEVM780103}{Normalized frequency of reverse turn, with weights (Levitt,
#'  1978)}
#'  \item{LEVM780104}{Normalized frequency of alpha-helix, unweighted (Levitt,
#'  1978)}
#'  \item{LEVM780105}{Normalized frequency of beta-sheet, unweighted (Levitt,
#'  1978)}
#'  \item{LEVM780106}{Normalized frequency of reverse turn, unweighted (Levitt,
#'  1978)}
#'  \item{LEWP710101}{Frequency of occurrence in beta-bends (Lewis et al., 1971)}
#'  \item{LIFS790101}{Conformational preference for all beta-strands
#'  (Lifson-Sander, 1979)}
#'  \item{LIFS790102}{Conformational preference for parallel beta-strands
#'  (Lifson-Sander, 1979)}
#'  \item{LIFS790103}{Conformational preference for antiparallel beta-strands
#'  (Lifson-Sander, 1979)}
#'  \item{MANP780101}{Average surrounding hydrophobicity (Manavalan-Ponnuswamy,
#'  1978)}
#'  \item{MAXF760101}{Normalized frequency of alpha-helix (Maxfield-Scheraga,
#'  1976)}
#'  \item{MAXF760102}{Normalized frequency of extended structure
#'  (Maxfield-Scheraga, 1976)}
#'  \item{MAXF760103}{Normalized frequency of zeta R (Maxfield-Scheraga, 1976)}
#'  \item{MAXF760104}{Normalized frequency of left-handed alpha-helix
#'  (Maxfield-Scheraga, 1976)}
#'  \item{MAXF760105}{Normalized frequency of zeta L (Maxfield-Scheraga, 1976)}
#'  \item{MAXF760106}{Normalized frequency of alpha region (Maxfield-Scheraga,
#'  1976)}
#'  \item{MCMT640101}{Refractivity (McMeekin et al., 1964), Cited by Jones (1975)}
#'  \item{MEEJ800101}{Retention coefficient in HPLC, pH7.4 (Meek, 1980)}
#'  \item{MEEJ800102}{Retention coefficient in HPLC, pH2.1 (Meek, 1980)}
#'  \item{MEEJ810101}{Retention coefficient in NaClO4 (Meek-Rossetti, 1981)}
#'  \item{MEEJ810102}{Retention coefficient in NaH2PO4 (Meek-Rossetti, 1981)}
#'  \item{MEIH800101}{Average reduced distance for C-alpha (Meirovitch et al.,
#'  1980)}
#'  \item{MEIH800102}{Average reduced distance for side chain (Meirovitch et al.,
#'  1980)}
#'  \item{MEIH800103}{Average side chain orientation angle (Meirovitch et al.,
#'  1980)}
#'  \item{MIYS850101}{Effective partition energy (Miyazawa-Jernigan, 1985)}
#'  \item{NAGK730101}{Normalized frequency of alpha-helix (Nagano, 1973)}
#'  \item{NAGK730102}{Normalized frequency of bata-structure (Nagano, 1973)}
#'  \item{NAGK730103}{Normalized frequency of coil (Nagano, 1973)}
#'  \item{NAKH900101}{AA composition of total proteins (Nakashima et al., 1990)}
#'  \item{NAKH900102}{SD of AA composition of total proteins (Nakashima et al.,
#'  1990)}
#'  \item{NAKH900103}{AA composition of mt-proteins (Nakashima et al., 1990)}
#'  \item{NAKH900104}{Normalized composition of mt-proteins (Nakashima et al.,
#'  1990)}
#'  \item{NAKH900105}{AA composition of mt-proteins from animal (Nakashima et al.,
#'  1990)}
#'  \item{NAKH900106}{Normalized composition from animal (Nakashima et al., 1990)}
#'  \item{NAKH900107}{AA composition of mt-proteins from fungi and plant (Nakashima
#'  et al., 1990)}
#'  \item{NAKH900108}{Normalized composition from fungi and plant (Nakashima et
#'  al., 1990)}
#'  \item{NAKH900109}{AA composition of membrane proteins (Nakashima et al., 1990)}
#'  \item{NAKH900110}{Normalized composition of membrane proteins (Nakashima et
#'  al., 1990)}
#'  \item{NAKH900111}{Transmembrane regions of non-mt-proteins (Nakashima et al.,
#'  1990)}
#'  \item{NAKH900112}{Transmembrane regions of mt-proteins (Nakashima et al.,
#'  1990)}
#'  \item{NAKH900113}{Ratio of average and computed composition (Nakashima et al.,
#'  1990)}
#'  \item{NAKH920101}{AA composition of CYT of single-spanning proteins
#'  (Nakashima-Nishikawa, 1992)}
#'  \item{NAKH920102}{AA composition of CYT2 of single-spanning proteins
#'  (Nakashima-Nishikawa, 1992)}
#'  \item{NAKH920103}{AA composition of EXT of single-spanning proteins
#'  (Nakashima-Nishikawa, 1992)}
#'  \item{NAKH920104}{AA composition of EXT2 of single-spanning proteins
#'  (Nakashima-Nishikawa, 1992)}
#'  \item{NAKH920105}{AA composition of MEM of single-spanning proteins
#'  (Nakashima-Nishikawa, 1992)}
#'  \item{NAKH920106}{AA composition of CYT of multi-spanning proteins
#'  (Nakashima-Nishikawa, 1992)}
#'  \item{NAKH920107}{AA composition of EXT of multi-spanning proteins
#'  (Nakashima-Nishikawa, 1992)}
#'  \item{NAKH920108}{AA composition of MEM of multi-spanning proteins
#'  (Nakashima-Nishikawa, 1992)}
#'  \item{NISK800101}{8 A contact number (Nishikawa-Ooi, 1980)}
#'  \item{NISK860101}{14 A contact number (Nishikawa-Ooi, 1986)}
#'  \item{NOZY710101}{Transfer energy, organic solvent/water (Nozaki-Tanford,
#'  1971)}
#'  \item{OOBM770101}{Average non-bonded energy per atom (Oobatake-Ooi, 1977)}
#'  \item{OOBM770102}{Short and medium range non-bonded energy per atom
#'  (Oobatake-Ooi, 1977)}
#'  \item{OOBM770103}{Long range non-bonded energy per atom (Oobatake-Ooi, 1977)}
#'  \item{OOBM770104}{Average non-bonded energy per residue (Oobatake-Ooi, 1977)}
#'  \item{OOBM770105}{Short and medium range non-bonded energy per residue
#'  (Oobatake-Ooi, 1977)}
#'  \item{OOBM850101}{Optimized beta-structure-coil equilibrium constant (Oobatake
#'  et al., 1985)}
#'  \item{OOBM850102}{Optimized propensity to form reverse turn (Oobatake et al.,
#'  1985)}
#'  \item{OOBM850103}{Optimized transfer energy parameter (Oobatake et al., 1985)}
#'  \item{OOBM850104}{Optimized average non-bonded energy per atom (Oobatake et
#'  al., 1985)}
#'  \item{OOBM850105}{Optimized side chain interaction parameter (Oobatake et al.,
#'  1985)}
#'  \item{PALJ810101}{Normalized frequency of alpha-helix from LG (Palau et al.,
#'  1981)}
#'  \item{PALJ810102}{Normalized frequency of alpha-helix from CF (Palau et al.,
#'  1981)}
#'  \item{PALJ810103}{Normalized frequency of beta-sheet from LG (Palau et al.,
#'  1981)}
#'  \item{PALJ810104}{Normalized frequency of beta-sheet from CF (Palau et al.,
#'  1981)}
#'  \item{PALJ810105}{Normalized frequency of turn from LG (Palau et al., 1981)}
#'  \item{PALJ810106}{Normalized frequency of turn from CF (Palau et al., 1981)}
#'  \item{PALJ810107}{Normalized frequency of alpha-helix in all-alpha class (Palau
#'  et al., 1981)}
#'  \item{PALJ810108}{Normalized frequency of alpha-helix in alpha+beta class
#'  (Palau et al., 1981)}
#'  \item{PALJ810109}{Normalized frequency of alpha-helix in alpha/beta class
#'  (Palau et al., 1981)}
#'  \item{PALJ810110}{Normalized frequency of beta-sheet in all-beta class (Palau
#'  et al., 1981)}
#'  \item{PALJ810111}{Normalized frequency of beta-sheet in alpha+beta class (Palau
#'  et al., 1981)}
#'  \item{PALJ810112}{Normalized frequency of beta-sheet in alpha/beta class (Palau
#'  et al., 1981)}
#'  \item{PALJ810113}{Normalized frequency of turn in all-alpha class (Palau et
#'  al., 1981)}
#'  \item{PALJ810114}{Normalized frequency of turn in all-beta class (Palau et al.,
#'  1981)}
#'  \item{PALJ810115}{Normalized frequency of turn in alpha+beta class (Palau et
#'  al., 1981)}
#'  \item{PALJ810116}{Normalized frequency of turn in alpha/beta class (Palau et
#'  al., 1981)}
#'  \item{PARJ860101}{HPLC parameter (Parker et al., 1986)}
#'  \item{PLIV810101}{Partition coefficient (Pliska et al., 1981)}
#'  \item{PONP800101}{Surrounding hydrophobicity in folded form (Ponnuswamy et al.,
#'  1980)}
#'  \item{PONP800102}{Average gain in surrounding hydrophobicity (Ponnuswamy et
#'  al., 1980)}
#'  \item{PONP800103}{Average gain ratio in surrounding hydrophobicity (Ponnuswamy
#'  et al., 1980)}
#'  \item{PONP800104}{Surrounding hydrophobicity in alpha-helix (Ponnuswamy et al.,
#'  1980)}
#'  \item{PONP800105}{Surrounding hydrophobicity in beta-sheet (Ponnuswamy et al.,
#'  1980)}
#'  \item{PONP800106}{Surrounding hydrophobicity in turn (Ponnuswamy et al., 1980)}
#'  \item{PONP800107}{Accessibility reduction ratio (Ponnuswamy et al., 1980)}
#'  \item{PONP800108}{Average number of surrounding residues (Ponnuswamy et al.,
#'  1980)}
#'  \item{PRAM820101}{Intercept in regression analysis (Prabhakaran-Ponnuswamy,
#'  1982)}
#'  \item{PRAM820102}{Slope in regression analysis x 1.0E1 (Prabhakaran-Ponnuswamy,
#'  1982)}
#'  \item{PRAM820103}{Correlation coefficient in regression analysis
#'  (Prabhakaran-Ponnuswamy, 1982)}
#'  \item{PRAM900101}{Hydrophobicity (Prabhakaran, 1990)}
#'  \item{PRAM900102}{Relative frequency in alpha-helix (Prabhakaran, 1990)}
#'  \item{PRAM900103}{Relative frequency in beta-sheet (Prabhakaran, 1990)}
#'  \item{PRAM900104}{Relative frequency in reverse-turn (Prabhakaran, 1990)}
#'  \item{PTIO830101}{Helix-coil equilibrium constant (Ptitsyn-Finkelstein, 1983)}
#'  \item{PTIO830102}{Beta-coil equilibrium constant (Ptitsyn-Finkelstein, 1983)}
#'  \item{QIAN880101}{Weights for alpha-helix at the window position of -6
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880102}{Weights for alpha-helix at the window position of -5
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880103}{Weights for alpha-helix at the window position of -4
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880104}{Weights for alpha-helix at the window position of -3
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880105}{Weights for alpha-helix at the window position of -2
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880106}{Weights for alpha-helix at the window position of -1
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880107}{Weights for alpha-helix at the window position of 0
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880108}{Weights for alpha-helix at the window position of 1
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880109}{Weights for alpha-helix at the window position of 2
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880110}{Weights for alpha-helix at the window position of 3
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880111}{Weights for alpha-helix at the window position of 4
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880112}{Weights for alpha-helix at the window position of 5
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880113}{Weights for alpha-helix at the window position of 6
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880114}{Weights for beta-sheet at the window position of -6
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880115}{Weights for beta-sheet at the window position of -5
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880116}{Weights for beta-sheet at the window position of -4
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880117}{Weights for beta-sheet at the window position of -3
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880118}{Weights for beta-sheet at the window position of -2
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880119}{Weights for beta-sheet at the window position of -1
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880120}{Weights for beta-sheet at the window position of 0
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880121}{Weights for beta-sheet at the window position of 1
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880122}{Weights for beta-sheet at the window position of 2
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880123}{Weights for beta-sheet at the window position of 3
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880124}{Weights for beta-sheet at the window position of 4
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880125}{Weights for beta-sheet at the window position of 5
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880126}{Weights for beta-sheet at the window position of 6
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880127}{Weights for coil at the window position of -6
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880128}{Weights for coil at the window position of -5
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880129}{Weights for coil at the window position of -4
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880130}{Weights for coil at the window position of -3
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880131}{Weights for coil at the window position of -2
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880132}{Weights for coil at the window position of -1
#'  (Qian-Sejnowski, 1988)}
#'  \item{QIAN880133}{Weights for coil at the window position of 0 (Qian-Sejnowski,
#'  1988)}
#'  \item{QIAN880134}{Weights for coil at the window position of 1 (Qian-Sejnowski,
#'  1988)}
#'  \item{QIAN880135}{Weights for coil at the window position of 2 (Qian-Sejnowski,
#'  1988)}
#'  \item{QIAN880136}{Weights for coil at the window position of 3 (Qian-Sejnowski,
#'  1988)}
#'  \item{QIAN880137}{Weights for coil at the window position of 4 (Qian-Sejnowski,
#'  1988)}
#'  \item{QIAN880138}{Weights for coil at the window position of 5 (Qian-Sejnowski,
#'  1988)}
#'  \item{QIAN880139}{Weights for coil at the window position of 6 (Qian-Sejnowski,
#'  1988)}
#'  \item{RACS770101}{Average reduced distance for C-alpha (Rackovsky-Scheraga,
#'  1977)}
#'  \item{RACS770102}{Average reduced distance for side chain (Rackovsky-Scheraga,
#'  1977)}
#'  \item{RACS770103}{Side chain orientational preference (Rackovsky-Scheraga,
#'  1977)}
#'  \item{RACS820101}{Average relative fractional occurrence in A0(i)
#'  (Rackovsky-Scheraga, 1982)}
#'  \item{RACS820102}{Average relative fractional occurrence in AR(i)
#'  (Rackovsky-Scheraga, 1982)}
#'  \item{RACS820103}{Average relative fractional occurrence in AL(i)
#'  (Rackovsky-Scheraga, 1982)}
#'  \item{RACS820104}{Average relative fractional occurrence in EL(i)
#'  (Rackovsky-Scheraga, 1982)}
#'  \item{RACS820105}{Average relative fractional occurrence in E0(i)
#'  (Rackovsky-Scheraga, 1982)}
#'  \item{RACS820106}{Average relative fractional occurrence in ER(i)
#'  (Rackovsky-Scheraga, 1982)}
#'  \item{RACS820107}{Average relative fractional occurrence in A0(i-1)
#'  (Rackovsky-Scheraga, 1982)}
#'  \item{RACS820108}{Average relative fractional occurrence in AR(i-1)
#'  (Rackovsky-Scheraga, 1982)}
#'  \item{RACS820109}{Average relative fractional occurrence in AL(i-1)
#'  (Rackovsky-Scheraga, 1982)}
#'  \item{RACS820110}{Average relative fractional occurrence in EL(i-1)
#'  (Rackovsky-Scheraga, 1982)}
#'  \item{RACS820111}{Average relative fractional occurrence in E0(i-1)
#'  (Rackovsky-Scheraga, 1982)}
#'  \item{RACS820112}{Average relative fractional occurrence in ER(i-1)
#'  (Rackovsky-Scheraga, 1982)}
#'  \item{RACS820113}{Value of theta(i) (Rackovsky-Scheraga, 1982)}
#'  \item{RACS820114}{Value of theta(i-1) (Rackovsky-Scheraga, 1982)}
#'  \item{RADA880101}{Transfer free energy from chx to wat (Radzicka-Wolfenden,
#'  1988)}
#'  \item{RADA880102}{Transfer free energy from oct to wat (Radzicka-Wolfenden,
#'  1988)}
#'  \item{RADA880103}{Transfer free energy from vap to chx (Radzicka-Wolfenden,
#'  1988)}
#'  \item{RADA880104}{Transfer free energy from chx to oct (Radzicka-Wolfenden,
#'  1988)}
#'  \item{RADA880105}{Transfer free energy from vap to oct (Radzicka-Wolfenden,
#'  1988)}
#'  \item{RADA880106}{Accessible surface area (Radzicka-Wolfenden, 1988)}
#'  \item{RADA880107}{Energy transfer from out to in(95\%buried)
#'  (Radzicka-Wolfenden, 1988)}
#'  \item{RADA880108}{Mean polarity (Radzicka-Wolfenden, 1988)}
#'  \item{RICJ880101}{Relative preference value at N" (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880102}{Relative preference value at N' (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880103}{Relative preference value at N-cap (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880104}{Relative preference value at N1 (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880105}{Relative preference value at N2 (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880106}{Relative preference value at N3 (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880107}{Relative preference value at N4 (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880108}{Relative preference value at N5 (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880109}{Relative preference value at Mid (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880110}{Relative preference value at C5 (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880111}{Relative preference value at C4 (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880112}{Relative preference value at C3 (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880113}{Relative preference value at C2 (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880114}{Relative preference value at C1 (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880115}{Relative preference value at C-cap (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880116}{Relative preference value at C' (Richardson-Richardson,
#'  1988)}
#'  \item{RICJ880117}{Relative preference value at C" (Richardson-Richardson,
#'  1988)}
#'  \item{ROBB760101}{Information measure for alpha-helix (Robson-Suzuki, 1976)}
#'  \item{ROBB760102}{Information measure for N-terminal helix (Robson-Suzuki,
#'  1976)}
#'  \item{ROBB760103}{Information measure for middle helix (Robson-Suzuki, 1976)}
#'  \item{ROBB760104}{Information measure for C-terminal helix (Robson-Suzuki,
#'  1976)}
#'  \item{ROBB760105}{Information measure for extended (Robson-Suzuki, 1976)}
#'  \item{ROBB760106}{Information measure for pleated-sheet (Robson-Suzuki, 1976)}
#'  \item{ROBB760107}{Information measure for extended without H-bond
#'  (Robson-Suzuki, 1976)}
#'  \item{ROBB760108}{Information measure for turn (Robson-Suzuki, 1976)}
#'  \item{ROBB760109}{Information measure for N-terminal turn (Robson-Suzuki,
#'  1976)}
#'  \item{ROBB760110}{Information measure for middle turn (Robson-Suzuki, 1976)}
#'  \item{ROBB760111}{Information measure for C-terminal turn (Robson-Suzuki,
#'  1976)}
#'  \item{ROBB760112}{Information measure for coil (Robson-Suzuki, 1976)}
#'  \item{ROBB760113}{Information measure for loop (Robson-Suzuki, 1976)}
#'  \item{ROBB790101}{Hydration free energy (Robson-Osguthorpe, 1979)}
#'  \item{ROSG850101}{Mean area buried on transfer (Rose et al., 1985)}
#'  \item{ROSG850102}{Mean fractional area loss (Rose et al., 1985)}
#'  \item{ROSM880101}{Side chain hydropathy, uncorrected for solvation (Roseman,
#'  1988)}
#'  \item{ROSM880102}{Side chain hydropathy, corrected for solvation (Roseman,
#'  1988)}
#'  \item{ROSM880103}{Loss of Side chain hydropathy by helix formation (Roseman,
#'  1988)}
#'  \item{SIMZ760101}{Transfer free energy (Simon, 1976), Cited by Charton-Charton
#'  (1982)}
#'  \item{SNEP660101}{Principal component I (Sneath, 1966)}
#'  \item{SNEP660102}{Principal component II (Sneath, 1966)}
#'  \item{SNEP660103}{Principal component III (Sneath, 1966)}
#'  \item{SNEP660104}{Principal component IV (Sneath, 1966)}
#'  \item{SUEM840101}{Zimm-Bragg parameter s at 20 C (Sueki et al., 1984)}
#'  \item{SUEM840102}{Zimm-Bragg parameter sigma x 1.0E4 (Sueki et al., 1984)}
#'  \item{SWER830101}{Optimal matching hydrophobicity (Sweet-Eisenberg, 1983)}
#'  \item{TANS770101}{Normalized frequency of alpha-helix (Tanaka-Scheraga, 1977)}
#'  \item{TANS770102}{Normalized frequency of isolated helix (Tanaka-Scheraga,
#'  1977)}
#'  \item{TANS770103}{Normalized frequency of extended structure (Tanaka-Scheraga,
#'  1977)}
#'  \item{TANS770104}{Normalized frequency of chain reversal R (Tanaka-Scheraga,
#'  1977)}
#'  \item{TANS770105}{Normalized frequency of chain reversal S (Tanaka-Scheraga,
#'  1977)}
#'  \item{TANS770106}{Normalized frequency of chain reversal D (Tanaka-Scheraga,
#'  1977)}
#'  \item{TANS770107}{Normalized frequency of left-handed helix (Tanaka-Scheraga,
#'  1977)}
#'  \item{TANS770108}{Normalized frequency of zeta R (Tanaka-Scheraga, 1977)}
#'  \item{TANS770109}{Normalized frequency of coil (Tanaka-Scheraga, 1977)}
#'  \item{TANS770110}{Normalized frequency of chain reversal (Tanaka-Scheraga,
#'  1977)}
#'  \item{VASM830101}{Relative population of conformational state A (Vasquez et
#'  al., 1983)}
#'  \item{VASM830102}{Relative population of conformational state C (Vasquez et
#'  al., 1983)}
#'  \item{VASM830103}{Relative population of conformational state E (Vasquez et
#'  al., 1983)}
#'  \item{VELV850101}{Electron-ion interaction potential (Veljkovic et al., 1985)}
#'  \item{VENT840101}{Bitterness (Venanzi, 1984)}
#'  \item{VHEG790101}{Transfer free energy to lipophilic phase (von
#'  Heijne-Blomberg, 1979)}
#'  \item{WARP780101}{Average interactions per side chain atom (Warme-Morgan,
#'  1978)}
#'  \item{WEBA780101}{RF value in high salt chromatography (Weber-Lacey, 1978)}
#'  \item{WERD780101}{Propensity to be buried inside (Wertz-Scheraga, 1978)}
#'  \item{WERD780102}{Free energy change of epsilon(i) to epsilon(ex)
#'  (Wertz-Scheraga, 1978)}
#'  \item{WERD780103}{Free energy change of alpha(Ri) to alpha(Rh) (Wertz-Scheraga,
#'  1978)}
#'  \item{WERD780104}{Free energy change of epsilon(i) to alpha(Rh)
#'  (Wertz-Scheraga, 1978)}
#'  \item{WOEC730101}{Polar requirement (Woese, 1973)}
#'  \item{WOLR810101}{Hydration potential (Wolfenden et al., 1981)}
#'  \item{WOLS870101}{Principal property value z1 (Wold et al., 1987)}
#'  \item{WOLS870102}{Principal property value z2 (Wold et al., 1987)}
#'  \item{WOLS870103}{Principal property value z3 (Wold et al., 1987)}
#'  \item{YUTK870101}{Unfolding Gibbs energy in water, pH7.0 (Yutani et al., 1987)}
#'  \item{YUTK870102}{Unfolding Gibbs energy in water, pH9.0 (Yutani et al., 1987)}
#'  \item{YUTK870103}{Activation Gibbs energy of unfolding, pH7.0 (Yutani et al.,
#'  1987)}
#'  \item{YUTK870104}{Activation Gibbs energy of unfolding, pH9.0 (Yutani et al.,
#'  1987)}
#'  \item{ZASB820101}{Dependence of partition coefficient on ionic strength
#'  (Zaslavsky et al., 1982)}
#'  \item{ZIMJ680101}{Hydrophobicity (Zimmerman et al., 1968)}
#'  \item{ZIMJ680102}{Bulkiness (Zimmerman et al., 1968)}
#'  \item{ZIMJ680103}{Polarity (Zimmerman et al., 1968)}
#'  \item{ZIMJ680104}{Isoelectric point (Zimmerman et al., 1968)}
#'  \item{ZIMJ680105}{RF rank (Zimmerman et al., 1968)}
#'  \item{AURR980101}{Normalized positional residue frequency at helix termini
#'  N4'(Aurora-Rose, 1998)}
#'  \item{AURR980102}{Normalized positional residue frequency at helix termini N"'
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980103}{Normalized positional residue frequency at helix termini N"
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980104}{Normalized positional residue frequency at helix termini
#'  N'(Aurora-Rose, 1998)}
#'  \item{AURR980105}{Normalized positional residue frequency at helix termini Nc
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980106}{Normalized positional residue frequency at helix termini N1
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980107}{Normalized positional residue frequency at helix termini N2
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980108}{Normalized positional residue frequency at helix termini N3
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980109}{Normalized positional residue frequency at helix termini N4
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980110}{Normalized positional residue frequency at helix termini N5
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980111}{Normalized positional residue frequency at helix termini C5
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980112}{Normalized positional residue frequency at helix termini C4
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980113}{Normalized positional residue frequency at helix termini C3
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980114}{Normalized positional residue frequency at helix termini C2
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980115}{Normalized positional residue frequency at helix termini C1
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980116}{Normalized positional residue frequency at helix termini Cc
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980117}{Normalized positional residue frequency at helix termini C'
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980118}{Normalized positional residue frequency at helix termini C"
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980119}{Normalized positional residue frequency at helix termini C"'
#'  (Aurora-Rose, 1998)}
#'  \item{AURR980120}{Normalized positional residue frequency at helix termini C4'
#'  (Aurora-Rose, 1998)}
#'  \item{ONEK900101}{Delta G values for the peptides extrapolated to 0 M urea
#'  (O'Neil-DeGrado, 1990)}
#'  \item{ONEK900102}{Helix formation parameters (delta delta G) (O'Neil-DeGrado,
#'  1990)}
#'  \item{VINM940101}{Normalized flexibility parameters (B-values), average
#'  (Vihinen et al., 1994)}
#'  \item{VINM940102}{Normalized flexibility parameters (B-values) for each residue
#'  surrounded by none rigid neighbours (Vihinen et al., 1994)}
#'  \item{VINM940103}{Normalized flexibility parameters (B-values) for each residue
#'  surrounded by one rigid neighbours (Vihinen et al., 1994)}
#'  \item{VINM940104}{Normalized flexibility parameters (B-values) for each residue
#'  surrounded by two rigid neighbours (Vihinen et al., 1994)}
#'  \item{MUNV940101}{Free energy in alpha-helical conformation (Munoz-Serrano,
#'  1994)}
#'  \item{MUNV940102}{Free energy in alpha-helical region (Munoz-Serrano, 1994)}
#'  \item{MUNV940103}{Free energy in beta-strand conformation (Munoz-Serrano,
#'  1994)}
#'  \item{MUNV940104}{Free energy in beta-strand region (Munoz-Serrano, 1994)}
#'  \item{MUNV940105}{Free energy in beta-strand region (Munoz-Serrano, 1994)}
#'  \item{WIMW960101}{Free energies of transfer of AcWl-X-LL peptides from bilayer
#'  interface to water (Wimley-White, 1996)}
#'  \item{KIMC930101}{Thermodynamic beta sheet propensity (Kim-Berg, 1993)}
#'  \item{MONM990101}{Turn propensity scale for transmembrane helices (Monne et
#'  al., 1999)}
#'  \item{BLAM930101}{Alpha helix propensity of position 44 in T4 lysozyme (Blaber
#'  et al., 1993)}
#'  \item{PARS000101}{p-Values of mesophilic proteins based on the distributions of
#'  B values (Parthasarathy-Murthy, 2000)}
#'  \item{PARS000102}{p-Values of thermophilic proteins based on the distributions
#'  of B values (Parthasarathy-Murthy, 2000)}
#'  \item{KUMS000101}{Distribution of amino acid residues in the 18 non-redundant
#'  families of thermophilic proteins (Kumar et al., 2000)}
#'  \item{KUMS000102}{Distribution of amino acid residues in the 18 non-redundant
#'  families of mesophilic proteins (Kumar et al., 2000)}
#'  \item{KUMS000103}{Distribution of amino acid residues in the alpha-helices in
#'  thermophilic proteins (Kumar et al., 2000)}
#'  \item{KUMS000104}{Distribution of amino acid residues in the alpha-helices in
#'  mesophilic proteins (Kumar et al., 2000)}
#'  \item{TAKK010101}{Side-chain contribution to protein stability (kJ/mol)
#'  (Takano-Yutani, 2001)}
#'  \item{FODM020101}{Propensity of amino acids within pi-helices
#'  (Fodje-Al-Karadaghi, 2002)}
#'  \item{NADH010101}{Hydropathy scale based on self-information values in the
#'  two-state model (5\% accessibility) (Naderi-Manesh et al., 2001)}
#'  \item{NADH010102}{Hydropathy scale based on self-information values in the
#'  two-state model (9\% accessibility) (Naderi-Manesh et al., 2001)}
#'  \item{NADH010103}{Hydropathy scale based on self-information values in the
#'  two-state model (16\% accessibility) (Naderi-Manesh et al., 2001)}
#'  \item{NADH010104}{Hydropathy scale based on self-information values in the
#'  two-state model (20\% accessibility) (Naderi-Manesh et al., 2001)}
#'  \item{NADH010105}{Hydropathy scale based on self-information values in the
#'  two-state model (25\% accessibility) (Naderi-Manesh et al., 2001)}
#'  \item{NADH010106}{Hydropathy scale based on self-information values in the
#'  two-state model (36\% accessibility) (Naderi-Manesh et al., 2001)}
#'  \item{NADH010107}{Hydropathy scale based on self-information values in the
#'  two-state model (50\% accessibility) (Naderi-Manesh et al., 2001)}
#'  \item{MONM990201}{Averaged turn propensities in a transmembrane helix (Monne et
#'  al., 1999)}
#'  \item{KOEP990101}{Alpha-helix propensity derived from designed sequences
#'  (Koehl-Levitt, 1999)}
#'  \item{KOEP990102}{Beta-sheet propensity derived from designed sequences
#'  (Koehl-Levitt, 1999)}
#'  \item{CEDJ970101}{Composition of amino acids in extracellular proteins
#'  (percent) (Cedano et al., 1997)}
#'  \item{CEDJ970102}{Composition of amino acids in anchored proteins (percent)
#'  (Cedano et al., 1997)}
#'  \item{CEDJ970103}{Composition of amino acids in membrane proteins (percent)
#'  (Cedano et al., 1997)}
#'  \item{CEDJ970104}{Composition of amino acids in intracellular proteins
#'  (percent) (Cedano et al., 1997)}
#'  \item{CEDJ970105}{Composition of amino acids in nuclear proteins (percent)
#'  (Cedano et al., 1997)}
#'  \item{FUKS010101}{Surface composition of amino acids in intracellular proteins
#'  of thermophiles (percent) (Fukuchi-Nishikawa, 2001)}
#'  \item{FUKS010102}{Surface composition of amino acids in intracellular proteins
#'  of mesophiles (percent) (Fukuchi-Nishikawa, 2001)}
#'  \item{FUKS010103}{Surface composition of amino acids in extracellular proteins
#'  of mesophiles (percent) (Fukuchi-Nishikawa, 2001)}
#'  \item{FUKS010104}{Surface composition of amino acids in nuclear proteins
#'  (percent) (Fukuchi-Nishikawa, 2001)}
#'  \item{FUKS010105}{Interior composition of amino acids in intracellular proteins
#'  of thermophiles (percent) (Fukuchi-Nishikawa, 2001)}
#'  \item{FUKS010106}{Interior composition of amino acids in intracellular proteins
#'  of mesophiles (percent) (Fukuchi-Nishikawa, 2001)}
#'  \item{FUKS010107}{Interior composition of amino acids in extracellular proteins
#'  of mesophiles (percent) (Fukuchi-Nishikawa, 2001)}
#'  \item{FUKS010108}{Interior composition of amino acids in nuclear proteins
#'  (percent) (Fukuchi-Nishikawa, 2001)}
#'  \item{FUKS010109}{Entire chain composition of amino acids in intracellular
#'  proteins of thermophiles (percent) (Fukuchi-Nishikawa, 2001)}
#'  \item{FUKS010110}{Entire chain composition of amino acids in intracellular
#'  proteins of mesophiles (percent) (Fukuchi-Nishikawa, 2001)}
#'  \item{FUKS010111}{Entire chain composition of amino acids in extracellular
#'  proteins of mesophiles (percent) (Fukuchi-Nishikawa, 2001)}
#'  \item{FUKS010112}{Entire chain compositino of amino acids in nuclear proteins
#'  (percent) (Fukuchi-Nishikawa, 2001)}
#'  \item{AVBF000101}{Screening coefficients gamma, local (Avbelj, 2000)}
#'  \item{AVBF000102}{Screening coefficients gamma, non-local (Avbelj, 2000)}
#'  \item{AVBF000103}{Slopes tripeptide, FDPB VFF neutral (Avbelj, 2000)}
#'  \item{AVBF000104}{Slopes tripeptides, LD VFF neutral (Avbelj, 2000)}
#'  \item{AVBF000105}{Slopes tripeptide, FDPB VFF noside (Avbelj, 2000)}
#'  \item{AVBF000106}{Slopes tripeptide FDPB VFF all (Avbelj, 2000)}
#'  \item{AVBF000107}{Slopes tripeptide FDPB PARSE neutral (Avbelj, 2000)}
#'  \item{AVBF000108}{Slopes dekapeptide, FDPB VFF neutral (Avbelj, 2000)}
#'  \item{AVBF000109}{Slopes proteins, FDPB VFF neutral (Avbelj, 2000)}
#'  \item{YANJ020101}{Side-chain conformation by gaussian evolutionary method (Yang
#'  et al., 2002)}
#'  \item{MITS020101}{Amphiphilicity index (Mitaku et al., 2002)}
#'  \item{TSAJ990101}{Volumes including the crystallographic waters using the
#'  ProtOr (Tsai et al., 1999)}
#'  \item{TSAJ990102}{Volumes not including the crystallographic waters using the
#'  ProtOr (Tsai et al., 1999)}
#'  \item{COSI940101}{Electron-ion interaction potential values (Cosic, 1994)}
#'  \item{PONP930101}{Hydrophobicity scales (Ponnuswamy, 1993)}
#'  \item{WILM950101}{Hydrophobicity coefficient in RP-HPLC, C18 with
#'  0.1\%TFA/MeCN/H2O (Wilce et al. 1995)}
#'  \item{WILM950102}{Hydrophobicity coefficient in RP-HPLC, C8 with
#'  0.1\%TFA/MeCN/H2O (Wilce et al. 1995)}
#'  \item{WILM950103}{Hydrophobicity coefficient in RP-HPLC, C4 with
#'  0.1\%TFA/MeCN/H2O (Wilce et al. 1995)}
#'  \item{WILM950104}{Hydrophobicity coefficient in RP-HPLC, C18 with
#'  0.1\%TFA/2-PrOH/MeCN/H2O (Wilce et al. 1995)}
#'  \item{KUHL950101}{Hydrophilicity scale (Kuhn et al., 1995)}
#'  \item{GUOD860101}{Retention coefficient at pH 2 (Guo et al., 1986)}
#'  \item{JURD980101}{Modified Kyte-Doolittle hydrophobicity scale (Juretic et al.,
#'  1998)}
#'  \item{BASU050101}{Interactivity scale obtained from the contact matrix
#'  (Bastolla et al., 2005)}
#'  \item{BASU050102}{Interactivity scale obtained by maximizing the mean of
#'  correlation coefficient over single-domain globular proteins (Bastolla et al.,
#'  2005)}
#'  \item{BASU050103}{Interactivity scale obtained by maximizing the mean of
#'  correlation coefficient over pairs of sequences sharing the TIM barrel fold
#'  (Bastolla et al., 2005)}
#'  \item{SUYM030101}{Linker propensity index (Suyama-Ohara, 2003)}
#'  \item{PUNT030101}{Knowledge-based membrane-propensity scale from 1D_Helix in
#'  MPtopo databases (Punta-Maritan, 2003)}
#'  \item{PUNT030102}{Knowledge-based membrane-propensity scale from 3D_Helix in
#'  MPtopo databases (Punta-Maritan, 2003)}
#'  \item{GEOR030101}{Linker propensity from all dataset (George-Heringa, 2003)}
#'  \item{GEOR030102}{Linker propensity from 1-linker dataset (George-Heringa,
#'  2003)}
#'  \item{GEOR030103}{Linker propensity from 2-linker dataset (George-Heringa,
#'  2003)}
#'  \item{GEOR030104}{Linker propensity from 3-linker dataset (George-Heringa,
#'  2003)}
#'  \item{GEOR030105}{Linker propensity from small dataset (linker length is less
#'  than six residues) (George-Heringa, 2003)}
#'  \item{GEOR030106}{Linker propensity from medium dataset (linker length is
#'  between six and 14 residues) (George-Heringa, 2003)}
#'  \item{GEOR030107}{Linker propensity from long dataset (linker length is greater
#'  than 14 residues) (George-Heringa, 2003)}
#'  \item{GEOR030108}{Linker propensity from helical (annotated by DSSP) dataset
#'  (George-Heringa, 2003)}
#'  \item{GEOR030109}{Linker propensity from non-helical (annotated by DSSP)
#'  dataset (George-Heringa, 2003)}
#'  \item{ZHOH040101}{The stability scale from the knowledge-based atom-atom
#'  potential (Zhou-Zhou, 2004)}
#'  \item{ZHOH040102}{The relative stability scale extracted from mutation
#'  experiments (Zhou-Zhou, 2004)}
#'  \item{ZHOH040103}{Buriability (Zhou-Zhou, 2004)}
#'  \item{BAEK050101}{Linker index (Bae et al., 2005)}
#'  \item{HARY940101}{Mean volumes of residues buried in protein interiors (Harpaz
#'  et al., 1994)}
#'  \item{PONJ960101}{Average volumes of residues (Pontius et al., 1996)}
#'  \item{DIGM050101}{Hydrostatic pressure asymmetry index, PAI (Di Giulio, 2005)}
#'  \item{WOLR790101}{Hydrophobicity index (Wolfenden et al., 1979)}
#'  \item{OLSK800101}{Average internal preferences (Olsen, 1980)}
#'  \item{KIDA850101}{Hydrophobicity-related index (Kidera et al., 1985)}
#'  \item{GUYH850102}{Apparent partition energies calculated from Wertz-Scheraga
#'  index (Guy, 1985)}
#'  \item{GUYH850103}{Apparent partition energies calculated from Robson-Osguthorpe
#'  index (Guy, 1985)}
#'  \item{GUYH850104}{Apparent partition energies calculated from Janin index (Guy,
#'  1985)}
#'  \item{GUYH850105}{Apparent partition energies calculated from Chothia index
#'  (Guy, 1985)}
#'  \item{ROSM880104}{Hydropathies of amino acid side chains, neutral form
#'  (Roseman, 1988)}
#'  \item{ROSM880105}{Hydropathies of amino acid side chains, pi-values in pH 7.0
#'  (Roseman, 1988)}
#'  \item{JACR890101}{Weights from the IFH scale (Jacobs-White, 1989)}
#'  \item{COWR900101}{Hydrophobicity index, 3.0 pH (Cowan-Whittaker, 1990)}
#'  \item{BLAS910101}{Scaled side chain hydrophobicity values (Black-Mould, 1991)}
#'  \item{CASG920101}{Hydrophobicity scale from native protein structures
#'  (Casari-Sippl, 1992)}
#'  \item{CORJ870101}{NNEIG index (Cornette et al., 1987)}
#'  \item{CORJ870102}{SWEIG index (Cornette et al., 1987)}
#'  \item{CORJ870103}{PRIFT index (Cornette et al., 1987)}
#'  \item{CORJ870104}{PRILS index (Cornette et al., 1987)}
#'  \item{CORJ870105}{ALTFT index (Cornette et al., 1987)}
#'  \item{CORJ870106}{ALTLS index (Cornette et al., 1987)}
#'  \item{CORJ870107}{TOTFT index (Cornette et al., 1987)}
#'  \item{CORJ870108}{TOTLS index (Cornette et al., 1987)}
#'  \item{MIYS990101}{Relative partition energies derived by the Bethe
#'  approximation (Miyazawa-Jernigan, 1999)}
#'  \item{MIYS990102}{Optimized relative partition energies - method A
#'  (Miyazawa-Jernigan, 1999)}
#'  \item{MIYS990103}{Optimized relative partition energies - method B
#'  (Miyazawa-Jernigan, 1999)}
#'  \item{MIYS990104}{Optimized relative partition energies - method C
#'  (Miyazawa-Jernigan, 1999)}
#'  \item{MIYS990105}{Optimized relative partition energies - method D
#'  (Miyazawa-Jernigan, 1999)}
#'  \item{ENGD860101}{Hydrophobicity index (Engelman et al., 1986)}
#'  \item{FASG890101}{Hydrophobicity index (Fasman, 1989)}
#'  \item{K6.5}{Values of Wc in proteins from class Beta, cutoff 6 A, 
#'  separation 5 (Wozniak, 2014)}
#'  \item{K8.5}{Values of Wc in proteins from class Beta, cutoff 8 A, 
#'  separation 5 (Wozniak, 2014)}
#'  \item{K12.5}{Values of Wc in proteins from class Beta, cutoff 12 A, 
#'  separation 5 (Wozniak, 2014)}
#'  \item{K6.15}{Values of Wc in proteins from class Beta, cutoff 6 A, 
#'  separation 15 (Wozniak, 2014)}
#'  \item{K8.15}{Values of Wc in proteins from class Beta, cutoff 8 A, 
#'  separation 15 (Wozniak, 2014)}
#'  \item{K12.15}{Values of Wc in proteins from class Beta, cutoff 12 A, 
#'  separation 15 (Wozniak, 2014)}
#' }
#' @references Kawashima, S. and Kanehisa, M. (2000) AAindex: amino acid 
#' index database. Nucleic Acids Res., 28:374.
#' 
#' Wozniak, P. and Kotulska M. (2014) Characteristics of protein 
#' residue-residue contacts and their application in contact prediction.
#' 20(11):2497
#' 
#' @source AAIndex database.
#' @keywords datasets
#' @examples 
#' data(aaprop)
#' 
NULL