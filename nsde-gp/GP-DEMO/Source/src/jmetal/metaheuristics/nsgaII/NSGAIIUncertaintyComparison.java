//  NSGAII.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.metaheuristics.nsgaII;

import java.io.File;
import java.util.Properties;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.Ranking;
import jmetal.util.RankingUnderUncertantyZaStPrimerjav;
import jmetal.core.*;
import jmetal.init.Params;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.*;

/**
 * This class implements the NSGA-II algorithm. 
 */
public class NSGAIIUncertaintyComparison extends Algorithm {

	
  Model modelGlobal_;
  int steviloVsehPrimerjav = 0;
  int[] steviloNOBBCorrect = new int[3];
  int[] steviloNOBBIncorrect = new int[6];
  int[] steviloCorrectBB = new int[3];
  int[] steviloIncorrectBB = new int[6];
  int steviloDodatnihIzracunov = 0;
  int steviloPrimerjavZDodatnimiIzracuni = 0;
  double sirinaIntervalaZaupanja =0;
  
  
  /**
   * Constructor
   * @param problem Problem to solve
   */
  public NSGAIIUncertaintyComparison(Problem problem) {
    super (problem) ;
  } // NSGAII
  
  /**
   * Constructor
   * @param problem Problem to solve
   * @param properties Properties of algorithm
   */
  public NSGAIIUncertaintyComparison(Problem problem, Properties properties, Model mod){
    super (problem, properties);  
    modelGlobal_ = mod;
  }

  /**   
   * Runs the NSGA-II algorithm.
   * @return a <code>SolutionSet</code> that is a set of non dominated solutions
   * as a result of the algorithm execution
   * @throws JMException 
   */
  public SolutionSet execute() throws JMException, ClassNotFoundException {
    int populationSize;
    int maxEvaluations;
    int evaluations;
    int frontGen;
    int frontFormat;
    
    inicializirajTebeloNaNic(steviloNOBBCorrect);
    inicializirajTebeloNaNic(steviloNOBBIncorrect);
    
    inicializirajTebeloNaNic(steviloCorrectBB);
    inicializirajTebeloNaNic(steviloIncorrectBB);
    

    String frontFileName, logFileName, genFileName, contents;
    QualityIndicator indicators;    // QualityIndicator object
    int requiredEvaluations;        // Use in the example of use of the
                                    // indicators object (see below)
    SolutionSet population;
    SolutionSet offspringPopulation;
    SolutionSet union;

    Operator mutationOperator;
    Operator crossoverOperator;
    Operator selectionOperator;

    Distance distance = new Distance();

    //Read the parameters
    populationSize = Integer.parseInt(getParameter("Algorithm.populationSize"));
    maxEvaluations = Integer.parseInt(getParameter("Algorithm.maxEvaluations"));
    
    indicators = (QualityIndicator) getInputParameter("indicators");
    
    frontGen=Integer.parseInt(getParameter("Output.frontGen"));
    frontFormat=Integer.parseInt(getParameter("Output.frontFormat"));
    
    frontFileName=getParameter("Output.frontFileName");
    logFileName = getParameter("Output.logFileName");
    genFileName = getParameter("Output.genFileName");
    
    //Initialize the variables
    // 1. Generiranje populacije staršev P(0)
    population = new SolutionSet(populationSize);   
    evaluations = 0;

    requiredEvaluations = 0;

    //Read the operators
    mutationOperator = operators_.get("mutation");
    crossoverOperator = operators_.get("crossover");
    selectionOperator = operators_.get("selection");

    
    Solution newSolution;
    
  //Naredimo nadomestni model  
    int stTockZaUcenje = Integer.parseInt(getParameter("Model.activeSet"));;
    SolutionSet setOfExactlyEvaluatedSolutions = new SolutionSet(stTockZaUcenje);    
    for (int i = 0; i < stTockZaUcenje; i++) {
      newSolution = new Solution(problem_);
      problem_.evaluate(newSolution);
      problem_.evaluateConstraints(newSolution);
      newSolution.setExactllyEvaluated(true);
      setOfExactlyEvaluatedSolutions.add(newSolution);
    }    
    modelGlobal_.update(setOfExactlyEvaluatedSolutions);
    
    
    //Params.setSeed(Params.getSeed() + 1);
    PseudoRandom.setRandomToNull();
    
 // Create the initial solutionSet
    // 1. Ovrednostimo začetno populacijo staršev
    for (int i = 0; i < populationSize; i++) {
      newSolution = new Solution(problem_);
      problem_.evaluate(newSolution);
      problem_.evaluateConstraints(newSolution);
      newSolution.setExactllyEvaluated(true);
      evaluations++;
      population.add(newSolution);      
    } //for       
    
    
    
    
    
    population.printVariablesAndObjectivesToFile(logFileName);   
    
    // Generations 
    while (evaluations < maxEvaluations) {  //Dokler ustavitveni kriterij ni izpolnjen
        
       System.out.println("Evaluations: "+evaluations); 
      // Create the offSpring solutionSet      
      // 2. Pripravimo prazno začetno populacijo potomcev Q(0)  
      offspringPopulation = new SolutionSet(populationSize);  
      
      // 4.8 Populacijo potomcev Q(t) generiraj iz populacije staršev P(t) z uporabo turnirske selekcije, križanja in mutacije
      Solution[] parents = new Solution[2];
      for (int i = 0; i < (populationSize / 2); i++) {
        if (evaluations < maxEvaluations) {
          //obtain parents
          parents[0] = (Solution) selectionOperator.execute(population);
          parents[1] = (Solution) selectionOperator.execute(population);
          Solution[] offSpring = (Solution[]) crossoverOperator.execute(parents);
          mutationOperator.execute(offSpring[0]);
          mutationOperator.execute(offSpring[1]);
          problem_.evaluate(offSpring[0]);
          problem_.evaluateConstraints(offSpring[0]);
          offSpring[0].setExactllyEvaluated(true);
          offSpring[0].setStandardDevianceToZero();
          problem_.evaluate(offSpring[1]);
          problem_.evaluateConstraints(offSpring[1]);
          offSpring[1].setExactllyEvaluated(true);
          offSpring[1].setStandardDevianceToZero();
          offspringPopulation.add(offSpring[0]);
          offspringPopulation.add(offSpring[1]);
          evaluations += 2;
        } // if                            
      } // for

      // Create the solutionSet union of solutionSet and offSpring
      // 4.2 Združi stari populaciji staršev in potomcev
      union = ((SolutionSet) population).union(offspringPopulation);

     
      SolutionSet kopijaUnion = kopirajDataSet(union);
      // Ranking the union 
      // 4.3 NE-DOMINIRANO SORTIRANJE
      
      //System.out.println("Pred Sta ok " + staSetaEnaka(kopijaUnion, union));
      
      //z BB
      RankingUnderUncertantyZaStPrimerjavZaBBox rankingZaBBox = new RankingUnderUncertantyZaStPrimerjavZaBBox(kopijaUnion, problem_, modelGlobal_, problem_ ); 
      steviloCorrectBB = sestejTabeli(steviloCorrectBB, rankingZaBBox.SteviloPravilnihPrimerjavZBB);
      steviloIncorrectBB = sestejTabeli(steviloIncorrectBB, rankingZaBBox.SteviloNepravilnihPrimerjavZBB);
      steviloDodatnihIzracunov = steviloDodatnihIzracunov + rankingZaBBox.steviloDodatnihEkzaktnihIzracunov;
      steviloPrimerjavZDodatnimiIzracuni = steviloPrimerjavZDodatnimiIzracuni + rankingZaBBox.steviloPrimerjavZDodatnimiIzracuni;
      sirinaIntervalaZaupanja =sirinaIntervalaZaupanja + rankingZaBBox.povprecnaSirinaIntervala;
      
      //kopijaUnion.get(0).setObjective(0, 12.34);
      //System.out.println("Sta seta ok pol " + staSetaEnaka(kopijaUnion, union));
      
      
    //brez BB       
      RankingUnderUncertantyZaStPrimerjav ranking  = new RankingUnderUncertantyZaStPrimerjav(union, problem_, modelGlobal_, problem_ );   //new Ranking(union);
      steviloVsehPrimerjav = steviloVsehPrimerjav + ranking.steviloPrimerjav;
      steviloNOBBCorrect = sestejTabeli(steviloNOBBCorrect, ranking.NOBBSteviloPravilnihPrimerjav);
      steviloNOBBIncorrect = sestejTabeli(steviloNOBBIncorrect, ranking.NOBBSteviloNepravilnihPrimerjav);
      
      /*int tempSt = 0;
      for(int i=0; i< ranking.NOBBSteviloNepravilnihPrimerjav.length; i++)
  			tempSt = tempSt + ranking.NOBBSteviloNepravilnihPrimerjav[i];
  			
      //System.out.println("Napake brezBB:  " + tempSt);
      System.out.print("Napake brezBB: ");
      for(int i=0; i< ranking.NOBBSteviloNepravilnihPrimerjav.length; i++)
      	System.out.print(ranking.NOBBSteviloNepravilnihPrimerjav[i] + "  ");
      System.out.println();
      
      //tempSt = 0;
      //for(int i=0; i< rankingZaBBox.SteviloNepravilnihPrimerjavZBB.length; i++)
  			//tempSt = tempSt + rankingZaBBox.SteviloNepravilnihPrimerjavZBB[i];
      //System.out.println("Napake zBB:  " + tempSt);
      
      System.out.print("Napake zBB: ");
      for(int i=0; i< rankingZaBBox.SteviloNepravilnihPrimerjavZBB.length; i++)
      	System.out.print(rankingZaBBox.SteviloNepravilnihPrimerjavZBB[i] + "  ");
      System.out.println();  */
     
      
      int remain = populationSize;
      int index = 0;
      SolutionSet front = null;
      // 4.4 Pripravimo novo prazno populacijo staršev 
      population.clear();

      // Obtain the next front
      front = ranking.getSubfront(index);
      
      //-------------------------------------ENVIRONMENTAL SELECTION--------------------------------------------
      // 4.5 V populacijo P daj prvih i front, ki še gredo cele v njo
      while ((remain > 0) && (remain >= front.size())) {
        //Assign crowding distance to individuals
        distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
        //Add the individuals of this front
        for (int k = 0; k < front.size(); k++) {
          population.add(front.get(k));
        } // for

        //Decrement remain
        remain = remain - front.size();

        //Obtain the next front
        index++;
        if (remain > 0) {
          front = ranking.getSubfront(index);
        } // if        
      } // while

      // Remain is less than front(index).size, insert only the best one
      // 4.6 Fronto F(i+1), ki ne gre več cela v populacijo P(t) sortiraj z uporabo metrike nakopičenosti
      // 4.7 Populacijo P(t) dopolni z osebki iz F(i+1), ki so najmanj nakopičeni
      
      if (remain > 0) {  // front contains individuals to insert                        
        distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
        front.sort(new CrowdingComparator());
        for (int k = 0; k < remain; k++) {
          population.add(front.get(k));
        } // for

        remain = 0;
      } // if                               
      //-------------------------------------ENVIRONMENTAL SELECTION--------------------------------------------
            
      // This piece of code shows how to use the indicator object into the code
      // of NSGA-II. In particular, it finds the number of evaluations required
      // by the algorithm to obtain a Pareto front with a hypervolume higher
      // than the hypervolume of the true Pareto front.
      if ((indicators != null) && (requiredEvaluations == 0)) {
        double HV = indicators.getHypervolume(population);
        if (HV >= (0.98 * indicators.getTrueParetoFrontHypervolume())) {
          requiredEvaluations = evaluations;
        } // if
      } // if
    
      /**
       * Write to the output file every frontGen generations
       */
      if(evaluations%frontGen == 0){ 
    	//Mihas
         ranking.getSubfront(0).printVariablesAndObjectivesToFile(frontFileName);
         //population.printGenerationsToFile(frontFileName, frontFormat);       
      }
      
      population.printVariablesAndObjectivesToFile(logFileName);  
      
      contents=String.format("%d\t\t", evaluations);
      population.printToFile(genFileName, contents);

      
      
    } // while

    
    //izpis statistike primerjav
    sirinaIntervalaZaupanja = sirinaIntervalaZaupanja/100;	//evaluations;
    //sirinaIntervalaZaupanja = sirinaIntervalaZaupanja/2;		//MODELIRAMO 2 * ZARADI SELEKCIJE
    
    String imeDatoteke = getParameter("Output.logFileName").replace("Log", "Primerjava");
    //imeDatoteke = "D:\\Doktorat\\Metka\\app\\jMetalCI\\Primerjave.txt";
    zapisHeaderjaZaRezultate(imeDatoteke);
    
    izpisRezultataPrimerjav(imeDatoteke, problem_.getName(), modelGlobal_.getName(), "brezBB", steviloVsehPrimerjav, steviloNOBBCorrect, steviloNOBBIncorrect, -1, -1, sirinaIntervalaZaupanja);    
    
    izpisRezultataPrimerjav(imeDatoteke, problem_.getName(), modelGlobal_.getName(), "BB", steviloVsehPrimerjav, steviloCorrectBB, steviloIncorrectBB, steviloDodatnihIzracunov, steviloPrimerjavZDodatnimiIzracuni, sirinaIntervalaZaupanja);
    
    
    try{
    	if(modelGlobal_.getName() != "RF")
    		modelGlobal_.closeRConnection();
    	
    	modelGlobal_.closeRConnection();
    }
    catch (Exception e) {
		}
    
    // Return as output parameter the required evaluations
    setOutputParameter("evaluations", requiredEvaluations);

    // Return the first non-dominated front
    // 4.9 Ovrednosti osebke iz populacije potomcev 
    Ranking ranking = new Ranking(population);
   
    return ranking.getSubfront(0);
  } // execute

  
	private String staSetaEnaka(SolutionSet kopija, SolutionSet union) {
	  
		for(int i = 0; i<kopija.size(); i++)
		{
			if(!staResitviEnaki(kopija.get(i), union.get(i)))
				return "ne";
		}
		
		return "da";
  }

	private boolean staResitviEnaki(Solution solution, Solution solution2) {
	  
		for(int i=0; i< solution.numberOfObjectives(); i++)
			if(solution.getObjective(i) != solution2.getObjective(i))
				return false;
		
	  return true;
  }

	private SolutionSet kopirajDataSet(SolutionSet union) {
	  
			SolutionSet ss = new SolutionSet(union.getMaxSize());
			
			for(int i=0; i< union.size();i++)
				ss.add(new Solution(union.get(i)));
			
			return ss;
  }

	private void zapisHeaderjaZaRezultate(String pot) throws JMException {
		SolutionSet s=new SolutionSet();
		
		File file = new File(pot);
				
		//System.out.println(file.length());
		
		if(file.length()==0)
				s.printToFile(pot, "Seed	Problem	Model	Relacije	skupaj	pravilno	napa�no	Pravilno:	p<q	q<p	p||q" + "	Napa�no:	q<p-namesto-p<q	p||q-namesto-p<q	" +
				"p<q-namesto-q<p	p||q-namesto-q<p	p<q-namesto-p||q	q<p-namesto-p||q	�t_dodatnih_primerjav	�t_dodatnih_izracunov	sirina_intervala");
	  
  }

	private void izpisRezultataPrimerjav(String pot, String imeProblema, String imeModela, String relacija, int stVsehPrimerjav, int[] stNOBBCorrect, int[] stNOBBIncorrect, int stDodatnihIzr, int stPrimerjavZDodatnimiIzracuni, double sirinaIntervalaZaupanja) throws JMException 
	{
		SolutionSet s=new SolutionSet();
		String nizPravilnihPrimerjav = "";
		String nizNapacnihPrimerjav = "";
		int stPravilnih = 0;
		int stNapacnih = 0;
		String seed;
		
		seed = String.valueOf(Params.getSeed());
		
		for(int i=0; i< stNOBBCorrect.length; i++)
		{
			nizPravilnihPrimerjav = nizPravilnihPrimerjav + stNOBBCorrect[i] + "	";
			stPravilnih = stPravilnih + stNOBBCorrect[i];
		}
		nizPravilnihPrimerjav = nizPravilnihPrimerjav.substring(0, nizPravilnihPrimerjav.length()-1);	
		
		for(int i=0; i< stNOBBIncorrect.length; i++)
		{
			nizNapacnihPrimerjav = nizNapacnihPrimerjav + stNOBBIncorrect[i] + "	";
			stNapacnih = stNapacnih + stNOBBIncorrect[i];
		}
		nizNapacnihPrimerjav = nizNapacnihPrimerjav.substring(0, nizNapacnihPrimerjav.length()-1);
		
		String stDodatnihPrimerjav = "/";		//pri BB imamo ta dodaten podatek
		if(stPrimerjavZDodatnimiIzracuni > 0)
			stDodatnihPrimerjav = String.valueOf(stPrimerjavZDodatnimiIzracuni);
		
		String stDodatnihIzracunov = "/";		//pri BB imamo ta dodaten podatek
		if(stDodatnihIzr > 0)
			stDodatnihIzracunov = String.valueOf(stDodatnihIzr);
		
		s.printToFile(pot,  seed + "	" + imeProblema + "	" + imeModela + "	" + relacija + "	" + String.valueOf(steviloVsehPrimerjav) + "	"+ stPravilnih + 
				"	"+ stNapacnih + "	|	" + nizPravilnihPrimerjav + "	|	" + nizNapacnihPrimerjav + "	" + stDodatnihPrimerjav + "	" + stDodatnihIzracunov + "	" + sirinaIntervalaZaupanja);
	  
  }

	private int[] sestejTabeli(int[] tabela1, int[] tabela2) {
	  
		for(int i=0; i< tabela1.length; i++)
			tabela1[i] = tabela1[i] + tabela2[i];
			
		return tabela1;
  }

	private void inicializirajTebeloNaNic(int[] tabela) {
	  
		for(int i=0; i<tabela.length; i++)
			tabela[i] = 0;
	  
  }
} // NSGA-II
