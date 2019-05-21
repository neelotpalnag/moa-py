//  Ranking.java
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

package jmetal.util;

import jmetal.core.Model;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.init.Params;
import jmetal.problems.BNH;
import jmetal.problems.ConstrEx;
import jmetal.problems.EKGSimulator;
import jmetal.problems.ExternEvaluator;
import jmetal.problems.Fonseca;
import jmetal.problems.Golinski;
import jmetal.problems.IntRealProblem;
import jmetal.problems.Kursawe;
import jmetal.problems.OKA1;
import jmetal.problems.OKA2;
import jmetal.problems.Osyczka;
import jmetal.problems.Osyczka2;
import jmetal.problems.Poloni;
import jmetal.problems.Schaffer;
import jmetal.problems.Srinivas;
import jmetal.problems.StoreSteel;
import jmetal.problems.Tanaka;
import jmetal.problems.Viennet2;
import jmetal.problems.Viennet3;
import jmetal.problems.Viennet4;
import jmetal.problems.Water;
import jmetal.problems.DTLZ.DTLZ1;
import jmetal.problems.DTLZ.DTLZ1a;
import jmetal.problems.DTLZ.DTLZ2;
import jmetal.problems.DTLZ.DTLZ3;
import jmetal.problems.DTLZ.DTLZ4;
import jmetal.problems.DTLZ.DTLZ5;
import jmetal.problems.DTLZ.DTLZ6;
import jmetal.problems.DTLZ.DTLZ7;
import jmetal.problems.LZ09.LZ09_F1;
import jmetal.problems.LZ09.LZ09_F2;
import jmetal.problems.LZ09.LZ09_F3;
import jmetal.problems.LZ09.LZ09_F4;
import jmetal.problems.LZ09.LZ09_F5;
import jmetal.problems.LZ09.LZ09_F6;
import jmetal.problems.LZ09.LZ09_F7;
import jmetal.problems.LZ09.LZ09_F8;
import jmetal.problems.LZ09.LZ09_F9;
import jmetal.problems.WFG.WFG1;
import jmetal.problems.WFG.WFG2;
import jmetal.problems.WFG.WFG3;
import jmetal.problems.WFG.WFG4;
import jmetal.problems.WFG.WFG5;
import jmetal.problems.WFG.WFG6;
import jmetal.problems.WFG.WFG7;
import jmetal.problems.WFG.WFG8;
import jmetal.problems.WFG.WFG9;
import jmetal.problems.ZDT.ZDT1;
import jmetal.problems.ZDT.ZDT2;
import jmetal.problems.ZDT.ZDT3;
import jmetal.problems.ZDT.ZDT4;
import jmetal.problems.ZDT.ZDT5;
import jmetal.problems.ZDT.ZDT6;
import jmetal.problems.cec2009Competition.CEC2009_UF1;
import jmetal.problems.cec2009Competition.CEC2009_UF10;
import jmetal.problems.cec2009Competition.CEC2009_UF2;
import jmetal.problems.cec2009Competition.CEC2009_UF3;
import jmetal.problems.cec2009Competition.CEC2009_UF4;
import jmetal.problems.cec2009Competition.CEC2009_UF5;
import jmetal.problems.cec2009Competition.CEC2009_UF6;
import jmetal.problems.cec2009Competition.CEC2009_UF7;
import jmetal.problems.cec2009Competition.CEC2009_UF8;
import jmetal.problems.cec2009Competition.CEC2009_UF9;
import jmetal.problems.singleObjective.Griewank;
import jmetal.problems.singleObjective.OneMax;
import jmetal.problems.singleObjective.Sphere;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.comparators.DominanceComparatorUnderUncertanty;
import jmetal.util.comparators.OverallConstraintViolationComparator;

import java.util.*;

/**
 * This class implements some facilities for ranking solutions.
 * Given a <code>SolutionSet</code> object, their solutions are ranked 
 * according to scheme proposed in NSGA-II; as a result, a set of subsets 
 * are obtained. The subsets are numbered starting from 0 (in NSGA-II, the 
 * numbering starts from 1); thus, subset 0 contains the non-dominated 
 * solutions, subset 1 contains the non-dominated solutions after removing those
 * belonging to subset 0, and so on.
 */
public class RankingUnderUncertantyZaStPrimerjav {
  
  /**
   * The <code>SolutionSet</code> to rank
   */
  private SolutionSet   solutionSet_ ;
  
  /**
   * An array containing all the fronts found during the search
   */
  private SolutionSet[] ranking_  ;
  
  /**
   * stores a <code>Comparator</code> for dominance checking
   */
  private static final Comparator dominance_ = new DominanceComparatorUnderUncertanty();
  
  /**
   * stores a <code>Comparator</code> for Overal Constraint Violation Comparator
   * checking
   */
  private static final Comparator constraint_ = new OverallConstraintViolationComparator();
  
  public int steviloPrimerjav;
  public int[] NOBBSteviloPravilnihPrimerjav = new int[3];
  public int[] NOBBSteviloNepravilnihPrimerjav = new int[6];
  Model modelGlobal_;
  static Problem problem;
  SolutionSet populationZaPrimerjanjePrimerjav;
  public double povprecnaSirinaIntervala;
    
  /** 
   * Constructor.
   * @param population The <code>SolutionSet</code> to be ranked.
   * If the solution with uncertainty is on the first front it has to be exactly evaluated
   * @throws JMException 
   * @throws ClassNotFoundException 
   */       
  public  RankingUnderUncertantyZaStPrimerjav(SolutionSet population, Problem problem_, Model model, Problem prob) throws JMException, ClassNotFoundException {        
    
  	steviloPrimerjav =0;
  	inicializirajTebeloNaNic(NOBBSteviloPravilnihPrimerjav);
  	inicializirajTebeloNaNic(NOBBSteviloNepravilnihPrimerjav);
  	povprecnaSirinaIntervala = 0;
  	
  	modelGlobal_ = model;
  	problem = createProblem(prob.getName());
  	populationZaPrimerjanjePrimerjav = new SolutionSet(200);
  	int neprimerljivi = 0;
  	
  	narediPrimernoKopijoPopulacijeZaNOBB(populationZaPrimerjanjePrimerjav, population);
  	
  	Solution s1;
  	Solution s2;

  	
  	
  	solutionSet_ = population ;
  	//int additionalEvaluationsNeeded = 0;

    // dominateMe[i] contains the number of solutions dominating i        
    int [] dominateMe = new int[solutionSet_.size()];

    // iDominate[k] contains the list of solutions dominated by k
    List<Integer> [] iDominate = new List[solutionSet_.size()];

    // front[i] contains the list of individuals belonging to the front i
    List<Integer> [] front = new List[solutionSet_.size()+1];
        
    // flagDominate is an auxiliar variable
    int flagDominate;
    
    int flagDominateNOBB;
    
    // couldBeDominated is an auxiliar variable
    boolean reevaluateSolutionAgain = false; 

    // Initialize the fronts 
    for (int i = 0; i < front.length; i++)
      front[i] = new LinkedList<Integer>();
        
    //-> Fast non dominated sorting algorithm
    for (int p = 0; p < solutionSet_.size(); p++) {
    	reevaluateSolutionAgain = false;
    // Initialice the list of individuals that i dominate and the number
    // of individuals that dominate me
      iDominate[p] = new LinkedList<Integer>();
      dominateMe[p] = 0;            
      // For all q individuals , calculate if p dominates q or vice versa
      for (int q = 0; q < solutionSet_.size(); q++) {
      if(p != q)
      	{
      	  steviloPrimerjav++;   
      	  //if(p==10 && q == 6)
      	  	//System.out.print("");
      	  
	        flagDominate =constraint_.compare(population.get(p),population.get(q));
	        if (flagDominate == 0) {
	          flagDominate =dominance_.compare(population.get(p),population.get(q));                                
	        }
	        
	        if (flagDominate == -1) {       //if P dominates Q
	          iDominate[p].add(new Integer(q)); //Add q to the set of solutions dominatet by p
	        } 
	        
	        else if (flagDominate == 1) {
	          dominateMe[p]++;  //Domination counter of p 
	        }
	        
	        else if(flagDominate == 0){
	        	neprimerljivi ++;
	        }
	        
	        else
	        {
	        	reevaluateSolutionAgain = true;
	        }

	        //za toènost primerjav
	        flagDominateNOBB =constraint_.compare(populationZaPrimerjanjePrimerjav.get(p),populationZaPrimerjanjePrimerjav.get(q));
	        if (flagDominateNOBB == 0) {
	          flagDominateNOBB =dominance_.compare(populationZaPrimerjanjePrimerjav.get(p),populationZaPrimerjanjePrimerjav.get(q));                                
	        }
	        
	        
	        if(flagDominate==flagDominateNOBB)
	        {
	        	if (flagDominate == -1) {       //if P dominates Q
	        		NOBBSteviloPravilnihPrimerjav[0]++;	//iDominate[p].add(new Integer(q)); //Add q to the set of solutions dominatet by p
		        } 
		        
		        else if (flagDominate == 1) {
		        	NOBBSteviloPravilnihPrimerjav[1]++;	//dominateMe[p]++;  //Domination counter of p 
		        }
		        
		        else if(flagDominate == 0){
		        	NOBBSteviloPravilnihPrimerjav[2]++;
		        }
	        	
	        }
	        
	        else 
	        {
		        	if (flagDominate == -1) {       //if P dominates Q
		        		if(flagDominateNOBB == 1)
		        			NOBBSteviloNepravilnihPrimerjav[0]++;		//napoved pove q<p v resnici pravilno je p<q
		        		else
		        			NOBBSteviloNepravilnihPrimerjav[1]++;
			        } 
			        
			        else if (flagDominate == 1) {
			        	if(flagDominateNOBB == -1)
		        			NOBBSteviloNepravilnihPrimerjav[2]++;		//napoved pove p<q v resnici pravilno je q<p
		        		else
		        			NOBBSteviloNepravilnihPrimerjav[3]++;
			        }
			        
			        else if(flagDominate == 0){
			        	if(flagDominateNOBB == -1)
			        	{
		        			NOBBSteviloNepravilnihPrimerjav[4]++;		//napoved pove p<q v resnici pravilno je p||q
		        			//System.out.println("NO: p = "+p+"  q="+q);
			        	}
		        		else
		        			NOBBSteviloNepravilnihPrimerjav[5]++;
			        }

	        }
	        
      	}      
      }
            
      if(dominateMe[p] == 0 && population.get(p).isExactllyEvaluated() == false && reevaluateSolutionAgain == true)
      {
      	problem_.evaluate(population.get(p)) ;
      	population.get(p).setStandardDevianceToZero();
        problem_.evaluateConstraints(population.get(p));        
        population.get(p).setExactllyEvaluated(true);        
      	p = -1;
      }
    }
    
    for(int p=0; p<solutionSet_.size(); p++) {
      // If nobody dominates p, p belongs to the first front
      if (dominateMe[p] == 0) {
        front[0].add(new Integer(p));
        population.get(p).setRank(0);
      }            
    }
        
    //Obtain the rest of fronts
    int i = 0;
    Iterator<Integer> it1, it2 ; // Iterators
    while (front[i].size()!= 0) {
      i++;
      it1 = front[i-1].iterator();
      while (it1.hasNext()) {
        it2 = iDominate[it1.next().intValue()].iterator();
        while (it2.hasNext()) {
          int index = it2.next().intValue();
          dominateMe[index]--;
          if (dominateMe[index]==0) {
            front[i].add(new Integer(index));
            solutionSet_.get(index).setRank(i);
          }
        }
      }
    }
    //<-
        
    ranking_ = new SolutionSet[i];
    //0,1,2,....,i-1 are front, then i fronts
    for (int j = 0; j < i; j++) {
      ranking_[j] = new SolutionSet(front[j].size());
      it1 = front[j].iterator();
      while (it1.hasNext()) {
                ranking_[j].add(population.get(it1.next().intValue()));       
      }
    }
    
    //return steviloPrimerjav;
   
  } // Ranking

  
  private void inicializirajTebeloNaNic(int[] tabela) {

  	for(int i=0; i<tabela.length; i++)
			tabela[i] = 0;
	  
  }


	private void narediPrimernoKopijoPopulacijeZaNOBB(SolutionSet populationZaPrimerjanje, SolutionSet population) throws ClassNotFoundException, JMException {

  	Solution s;

  	/*for(int i=0; i< population.size()/2;i++)
  	{
  		s = new Solution(population.get(i));
  		populationZaPrimerjanje.add(s);
  	}*/
  	
  	for(int i=0; i< population.size();i++)
  	{
  		s = new Solution(population.get(i));
  		modelGlobal_.evaluate(s);
    	problem.evaluateConstraints(s);
  		s.setStandardDevianceToZero();
    	s.setExactllyEvaluated(true);
  		populationZaPrimerjanje.add(s);
  	}
  }


	public static Problem createProblem(String problemName) throws ClassNotFoundException{
      
  	int numberOfVariables = 2;
  	
      if (problemName.compareTo("DTLZ1") == 0) {          // 
          problem = new DTLZ1("Real");
      } else if (problemName.compareTo("DTLZ1a") == 0) {          // 
          problem = new DTLZ1a("Real");    
      } else if (problemName.compareTo("DTLZ2") == 0) {    // DTLZ2
          problem = new DTLZ2("Real");
      } else if (problemName.compareTo("DTLZ3") == 0) {    // DTLZ3  
          problem = new DTLZ3("Real");
      } else if (problemName.compareTo("DTLZ4") == 0) {    // DTLZ4
          problem = new DTLZ4("Real");
      } else if (problemName.compareTo("DTLZ5") == 0) {    // DTLZ5
          problem = new DTLZ5("Real");
      } else if (problemName.compareTo("DTLZ6") == 0) {    // DTLZ6
          problem = new DTLZ6("Real");
      } else if (problemName.compareTo("DTLZ7") == 0) {    // DTLZ7
          problem = new DTLZ7("Real", 23, 4);
      } else if (problemName.compareTo("LZ09") == 0) {    // LZ09
          // ni podrazred razreda Problem
          // LZ09 (int nvar, int nobj, int ptype, int dtype, int ltype)
      } else if (problemName.compareTo("LZ09_F1") == 0) {   // LZ09_F1
          problem = new LZ09_F1("Real");
      } else if (problemName.compareTo("LZ09_F2") == 0) {  // LZ09_F2
          problem = new LZ09_F2("Real");
      } else if (problemName.compareTo("LZ09_F3") == 0) {  // LZ09_F3
          problem = new LZ09_F3("Real");
      } else if (problemName.compareTo("LZ09_F4") == 0) {  // LZ09_F4
          problem = new LZ09_F4("Real");
      } else if (problemName.compareTo("LZ09_F5") == 0) {  // LZ09_F5
          problem = new LZ09_F5("Real");
      } else if (problemName.compareTo("LZ09_F6") == 0) {  // LZ09_F6
          problem = new LZ09_F6("Real");
      } else if (problemName.compareTo("LZ09_F7") == 0) {  // LZ09_F7
          problem = new LZ09_F7("Real");
      } else if (problemName.compareTo("LZ09_F8") == 0) {  // LZ09_F8
          problem = new LZ09_F8("Real");
      } else if (problemName.compareTo("LZ09_F9") == 0) {  // LZ09_F9
          problem = new LZ09_F9("Real");
      } else if (problemName.compareTo("Shapes") == 0) {   // Shapes
          // ni podrazred razreda Problem
      } else if (problemName.compareTo("Transformations") == 0) {  // Transformations
          // ni podrazred razreda Problem
      } else if (problemName.compareTo("WFG") == 0) {  // WFG
          // abstract class
      } else if (problemName.compareTo("WFG1") == 0) { // WFG1
          problem = new WFG1("Real");												
      } else if (problemName.compareTo("WFG2") == 0) { // WFG2
          problem = new WFG2("Real");
      } else if (problemName.compareTo("WFG3") == 0) { // WFG3
          problem = new WFG3("Real");
      } else if (problemName.compareTo("WFG4") == 0) { // WFG4
          problem = new WFG4("Real");
      } else if (problemName.compareTo("WFG5") == 0) { // WFG5
          problem = new WFG5("Real");
      } else if (problemName.compareTo("WFG6") == 0) { // WFG6
          problem = new WFG6("Real");
      } else if (problemName.compareTo("WFG7") == 0) { // WFG7
          problem = new WFG7("Real");
      } else if (problemName.compareTo("WFG8") == 0) { // WFG8
          problem = new WFG8("Real");
      } else if (problemName.compareTo("WFG9") == 0) { // WFG9
          problem = new WFG9("Real");
      } else if (problemName.compareTo("ZDT1") == 0) { // ZDT1
          problem = new ZDT1("Real");
      } else if (problemName.compareTo("ZDT2") == 0) { // ZDT2
          problem = new ZDT2("Real");
      } else if (problemName.compareTo("ZDT3") == 0) { // ZDT3
          problem = new ZDT3("Real");
      } else if (problemName.compareTo("ZDT4") == 0) { // ZDT4
          problem = new ZDT4("Real");
      } else if (problemName.compareTo("ZDT5") == 0) { // ZDT5
          problem = new ZDT5("Real");
      } else if (problemName.compareTo("ZDT6") == 0) { // ZDT6
          problem = new ZDT6("Real");
      } else if (problemName.compareTo("CEC2009_UF1") == 0) { // CEC2009_UF1
          problem = new CEC2009_UF1("Real");
      } else if (problemName.compareTo("CEC2009_UF2") == 0) { // CEC2009_UF2
          problem = new CEC2009_UF2("Real");
      } else if (problemName.compareTo("CEC2009_UF3") == 0) { // CEC2009_UF3
          problem = new CEC2009_UF3("Real");
      } else if (problemName.compareTo("CEC2009_UF4") == 0) { // CEC2009_UF4
          problem = new CEC2009_UF4("Real");
      } else if (problemName.compareTo("CEC2009_UF5") == 0) { // CEC2009_UF5
          problem = new CEC2009_UF5("Real");
      } else if (problemName.compareTo("CEC2009_UF6") == 0) { // CEC2009_UF6
          problem = new CEC2009_UF6("Real");
      } else if (problemName.compareTo("CEC2009_UF7") == 0) { // CEC2009_UF7
          problem = new CEC2009_UF7("Real");
      } else if (problemName.compareTo("CEC2009_UF8") == 0) { // CEC2009_UF8
          problem = new CEC2009_UF8("Real");
      } else if (problemName.compareTo("CEC2009_UF9") == 0) { // CEC2009_UF9
          problem = new CEC2009_UF9("Real");
      } else if (problemName.compareTo("CEC2009_UF10") == 0) { // CEC2009_UF10
          problem = new CEC2009_UF10("Real");
      } else if (problemName.compareTo("Griewank") == 0) { // Griewank
          problem = new Griewank("Real", numberOfVariables);
      } else if (problemName.compareTo("OneMax") == 0) { // OneMax
          // OneMax(Integer numberOfBits)
          problem = new OneMax(1);
      } else if (problemName.compareTo("Sphere") == 0) { // Sphere
          problem = new Sphere("Real", numberOfVariables);
      } else if (problemName.compareTo("TSP") == 0) { // TSP
          // TSP(String filename) throws FileNotFoundException, IOException, ClassNotFoundException
          // It accepts data files from TSPLIB (The file containing the definition of the problem)
          // problem = new TSP("");
      } else if (problemName.compareTo("ConstrEx") == 0) { // ConstrEx
          problem = new ConstrEx("Real");
      } else if (problemName.compareTo("Fonseca") == 0) { // Fonseca
          problem = new Fonseca("Real");
      } else if (problemName.compareTo("Golinski") == 0) { // Golinski
          problem = new Golinski("Real");
      } else if (problemName.compareTo("IntRealProblem") == 0) { // IntRealProblem
          problem = new IntRealProblem("IntReal");
      } else if (problemName.compareTo("Kursawe") == 0) { // Kursawe
          problem = new Kursawe("Real", numberOfVariables);
      } else if (problemName.compareTo("OKA1") == 0) { // OKA1
          problem = new OKA1("Real");
      } else if (problemName.compareTo("OKA2") == 0) { // OKA2
          problem = new OKA2("Real");
      } else if (problemName.compareTo("Osyczka") == 0) { // Osyczka2
          problem = new Osyczka("Real");
      } else if (problemName.compareTo("Osyczka2") == 0) { // Osyczka2
        problem = new Osyczka2("Real");
      } else if (problemName.compareTo("Poloni") == 0) { // Poloni
          problem = new Poloni("Real");
      } else if (problemName.compareTo("ProblemFactory") == 0) { // ProblemFactory
          // Class ProblemFactory represents a factory for problems, not a problem.
      } else if (problemName.compareTo("Schaffer") == 0) { // Schaffer
          problem = new Schaffer("Real");
      } else if (problemName.compareTo("Srinivas") == 0) { // Srinivas
          problem = new Srinivas("Real");
      } else if (problemName.compareTo("Tanaka") == 0) { // Tanaka
          problem = new Tanaka("Real");
      } else if (problemName.compareTo("Viennet2") == 0) { // Viennet2
          problem = new Viennet2("Real");
      } else if (problemName.compareTo("Viennet3") == 0) { // Viennet3
          problem = new Viennet3("Real");
      } else if (problemName.compareTo("Viennet4") == 0) { // Viennet4
          problem = new Viennet4("Real");
      } else if (problemName.compareTo("Water") == 0) { // Water
          problem = new Water("Real");
      } else if (problemName.compareTo("BNH") == 0) { // Water
        problem = new BNH("Real");
      } else if (problemName.compareTo("Extern_Evaluator") == 0) { // Extern_Evaluator	
          problem = new ExternEvaluator();
      } else if (problemName.compareTo("StoreSteel") == 0) { // Extern_Evaluator	
        problem = new StoreSteel();
        double [] topBorders = {1, 2, 15};
        Params.setFrontForHypervolume(topBorders);
      } else if (problemName.compareTo("EKGSimulator") == 0) { // Extern_Evaluator	
        problem = new EKGSimulator();
        double [] topBorders = {1, 2, 15};
        Params.setFrontForHypervolume(topBorders);
      } else {
          System.out.printf("Name %s of the PROBLEM in INPUT.txt is invalid", problemName);
      }

      double [] topBorders = new double[problem.getNumberOfObjectives()];
      for(int i=0; i<problem.getNumberOfObjectives(); i++)
      	topBorders[i] = 1000;
      
      if(problem.getName() == "StoreSteel")
      	topBorders = new double[]{1, 2, 15};
      
      else if(problem.getName() == "ExternEvaluator")
      	topBorders = Params.getLimitsForHypervolume();
      
      else if(problem.getName() == "Sphere")
      	topBorders = new double[]{2};
      	
      else if(problem.getName() == "Fonseca")
      	topBorders = new double[]{1, 1};
      
      else if(problem.getName() == "Golinski")
      	topBorders = new double[]{6000, 2000};
      
      else if(problem.getName() == "Kursawe")
      	topBorders = new double[]{-31, 2}; 	//-21, 2
      
      else if(problem.getName() == "Osyczka")
      	topBorders = new double[]{0, 80};
      
      else if(problem.getName() == "Osyczka2")
      	topBorders = new double[]{0, 80};
      
      else if(problem.getName() == "Poloni")
      	topBorders = new double[]{18, 30};
      
      else if(problem.getName() == "Srinivas")
      	topBorders = new double[]{250, 50};
      
      else if(problem.getName() == "Tanaka")
      	topBorders = new double[]{2, 2};
      
      else if(problem.getName() == "Viennet2")
      	topBorders = new double[]{5, -18, -14};
      	
      else if(problem.getName().contains("WFG"))
      {
      	for(int i=0; i<problem.getNumberOfObjectives(); i++)
        	topBorders[i] = 10;
      }
      
      else if(problem.getName().contains("DTLZ3"))
      	topBorders = new double[]{100, 100, 100};
      
      else if(problem.getName().contains("DTLZ"))
      {
      	for(int i=0; i<problem.getNumberOfObjectives(); i++)
        	topBorders[i] = 10;
      }
      
      else if(problem.getName().contains("BNH"))
      	topBorders = new double[]{150, 50};
      
      else if(problem.getName().contains("EKGSimulator"))
      	topBorders = new double[]{2, 2};
      
      
      Params.setFrontForHypervolume(topBorders);
      
      //if (problemName.compareTo("Extern_Evaluator") != 0)
      //	problem.setNumberOfFeatures(problem.getNumberOfObjectives());
      
      return problem;
  }

  
	private String[] approximateSolution(Solution child) throws ClassNotFoundException, JMException {

  	double distanceToFisibility;
  	double standardDevianceOfGlobalModel;
  	double[] confidenceIntervalsOfGlobalModel = new double[child.numberOfObjectives()];
  	double[] confidenceIntervalsOfLocalModel = new double[child.numberOfObjectives()];
  	double[] approximatedObjectivesOfGlobalModel = new double[child.numberOfObjectives()];
  	double[] approximatedObjectivesOfLocalModel = new double[child.numberOfObjectives()];
  	String[] allObjectivesAndConfidences = new String[child.numberOfObjectives() * 2 * 2  + 4];	//obj + confidence, 2 models, 4 labels
  	
  	Solution globalModelSolution;
  	
  	
  	
  	modelGlobal_.evaluate(child);
  	problem.evaluateConstraints(child);
  	//write data for log file
  	for(int i=0; i<child.numberOfObjectives(); i++)
		{
  		confidenceIntervalsOfGlobalModel[i] = child.getstandardDeviance(i);
  		approximatedObjectivesOfGlobalModel[i] = child.getObjective(i);
		}
  	
  	child.setExactllyEvaluated(false);
  	
  	
  		return allObjectivesAndConfidences;
  }

	
  /**
   * Returns a <code>SolutionSet</code> containing the solutions of a given rank. 
   * @param rank The rank
   * @return Object representing the <code>SolutionSet</code>.
   */
  public SolutionSet getSubfront(int rank) {
    return ranking_[rank];
  } // getSubFront

  /** 
  * Returns the total number of subFronts founds.
  */
  public int getNumberOfSubfronts() {
    return ranking_.length;
  } // getNumberOfSubfronts
} // Ranking
