//  SPGP.java
//
//  Author:
//       Dejan Petelin <dejan.petelin@ijs.si>
//
//  Copyright (c) 2011 Dejan Petelin
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

package jmetal.models;

import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import jmetal.core.Model;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.core.Variable;
import jmetal.util.JMException;
import spgp.SPGP;
import sun.misc.BASE64Encoder;

import com.mathworks.toolbox.javabuilder.MWCharArray;
import com.mathworks.toolbox.javabuilder.MWClassID;
import com.mathworks.toolbox.javabuilder.MWException;
import com.mathworks.toolbox.javabuilder.MWNumericArray;
import com.sun.java.swing.plaf.windows.resources.windows;

/**
 * Class representing model SPGP (Sparse Gaussian Processes using pseudo-inputs)
 */
public class SPGPm extends Model {
  
  private static final long serialVersionUID = 6970435628609740586L;
	private SPGP spgp_;
	private boolean fixed_;
	private int hypLength_;
	private MWNumericArray[] hyp_;
	private MWNumericArray[] bvs_;
	private ArrayList<double[]>[] inputs_;
	private ArrayList<Double>[] outputs_;
	private int numberOfPseudoInputs_;
	private List<String> data_;
	private boolean resetHyp_;
	private int maxOptIter_;

	private int windowSize_;
	private boolean resetWin_;
	private ArrayList<double[]> inputsAll_;
	private ArrayList<Double>[] outputsAll_;
	
	private ArrayList<Double>[] pseudoOutputs_;

  /** 
   * Constructor.
   * Creates a default instance of the SPGP model.
   * @throws Exception 
   * @throws ClassNotFoundException
   */
  public SPGPm(int noInput, int noOutput, int noPseudo) throws Exception {
  	super(noInput, noOutput);
  	
  	double[][] hyp;
  	
  	modelName_ = "SPGP";
  	resetHyp_ = true;
  	maxOptIter_ = -500;
  	numberOfPseudoInputs_ = noPseudo;
  	windowSize_ = noPseudo;
  	resetWin_ = false;
  	// SPGP models
  	spgp_ = new SPGP();
  	spgp_.init(); // random seed set
  	hyp_ = new MWNumericArray[noOutput];
  	bvs_ = new MWNumericArray[noOutput];
  	// get number of hyperparameter according to covariance function and number of inputs
  	hypLength_ = ((MWNumericArray)spgp_.hypLen(1, noInput)[0]).getInt();
  	// set hyperparameter values
		hyp = new double[hypLength_][];
  	for (int i = 0; i < hypLength_; i++)
		  hyp[i] = new double[]{0};
		for (int i = 0; i < noOutput; i++) {
  		hyp_[i] = new MWNumericArray(hyp); // default hyperparameter values
  		bvs_[i] = new MWNumericArray(MWClassID.DOUBLE);
		}
		fixed_ = false;
  	// set data
  	inputs_ = (ArrayList<double[]>[])new ArrayList[numberOfOutputs_];
  	outputs_ = (ArrayList<Double>[])new ArrayList[numberOfOutputs_];
  	outputsAll_ = (ArrayList<Double>[])new ArrayList[numberOfOutputs_];
  	pseudoOutputs_ = (ArrayList<Double>[])new ArrayList[numberOfOutputs_];
  	for (int i = 0; i < numberOfOutputs_; i++) {
  	  inputs_[i] = new ArrayList<double[]>();
  	  outputs_[i] = new ArrayList<Double>();
  	  outputsAll_[i] = new ArrayList<Double>();
  	  pseudoOutputs_[i] = new ArrayList<Double>();
  	}
	  inputsAll_ = new ArrayList<double[]>();
  	data_ = new ArrayList<String>();
  } // SPGP

  public void setResetHyp(boolean resetHyp) {
  	// whether or not reset hyperparameter values
  	resetHyp_ = resetHyp;
  }

  public void setWindowSize(int windowSize) {
  	if (resetWin_ && windowSize < numberOfPseudoInputs_)
  		windowSize_ = numberOfPseudoInputs_;
  	else
  		windowSize_ = windowSize;
  }

  /*public void setResetWin(boolean resetWin) {
  	// whether or not include BV in training window
  	resetWin_ = resetWin;
  	setWindowSize(windowSize_);
  }*/

  public void setHyperparameters(double[][] hyp) throws Exception {
		if (hyp.length != hypLength_ || hyp[0].length != numberOfOutputs_)
			throw new Exception("Dimenzija podanih hiperparametrov ni OK");
		double[][] hypt;
		for (int i = 0; i < numberOfOutputs_; i++) {
			hypt = new double[hypLength_][];
			for (int j = 0; j < hypLength_; j++)
				hypt[j] = new double[]{hyp[j][i]};
  		hyp_[i] = new MWNumericArray(hypt); // default hyperparameter values
  		bvs_[i] = new MWNumericArray(MWClassID.DOUBLE);
		}
		fixed_ = true;  	
  }

  public void setMaxOptIter(int maxOptIter) {
  	maxOptIter_ = maxOptIter;
  }

	@Override
  public void update(SolutionSet solutions) throws JMException, ClassNotFoundException {
		update(solutions, false);
	}
	
  public void update(SolutionSet solutions, boolean reset) throws JMException, ClassNotFoundException {
		Solution solution;
		
		// clear current inputs and outputs
		for (int i = 0; i < numberOfOutputs_; i++) {
		  inputs_[i].clear();
		  outputs_[i].clear();
		}
		for (int i = 0; i < solutions.size(); i++) {
			solution = solutions.get(i);
			append(solution);
		}
		try {
			optimize(reset);
    }
		catch (Exception e) {
	    e.printStackTrace();
    }
  }

	@Override
  public void update(Solution solution) throws JMException, ClassNotFoundException {
		update(solution, false);
	}
	
  public void update(Solution solution, boolean reset) throws JMException, ClassNotFoundException {
		// clear current inputs and outputs
		for (int i = 0; i < numberOfOutputs_; i++) {
		  inputs_[i].clear();
		  outputs_[i].clear();
		}
		append(solution);
		try {
	    optimize(reset);
    }
		catch (Exception e) {
	    e.printStackTrace();
    }
  }

	@Override
  public void evaluate(SolutionSet solutions) throws JMException, ClassNotFoundException {
		for (int i = 0; i < solutions.size(); i++) {
			evaluate(solutions.get(i));
		}
  }

	@Override
  public void evaluate(Solution solution) throws JMException, ClassNotFoundException {
		Object[] output;
		double variance;
		
		for (int i = 0; i < numberOfOutputs_; i++) {
			try {
				//output = spgp_.pred(3, toArray(outputs_[i]), toArray(inputs_[i]), bvs_[i], toArray(solution.getDecisionVariables()), hyp_[i], 0, 0);
				output = spgp_.pred(3, toArray(outputsAll_[i]), toArray(inputsAll_), bvs_[i], toArray(solution.getDecisionVariables()), hyp_[i], 0, 0);
				if (((MWCharArray)output[2]).toString().equalsIgnoreCase("OK")) {
		      if (solution.numberOfFeatures() > 0)
		      	solution.setFeature(i, ((MWNumericArray)output[0]).getDouble());
		      else
		      	solution.setObjective(i, ((MWNumericArray)output[0]).getDouble());
		      variance = ((MWNumericArray)output[1]).getDouble();
		      solution.setStandardDeviance(i, Math.sqrt(variance) * 3.0);
				}
				else {
		      if (solution.numberOfFeatures() > 0)
		      	solution.setFeature(i, Double.NaN);
		      else
		      	solution.setObjective(i, Double.NaN);
		      solution.setStandardDeviance(i, Double.NaN);
		      //e.printStackTrace();
		      System.out.println("Napaka pri evaluaciji!");
	      }
			}
			catch (MWException e) {
				System.out.println("[!] Napaka pri evaluaciji: " + e.getMessage() + "!");
			}
      solution.setExactllyEvaluated(false);
		}
  }

	private void append(Solution solution) throws JMException {
		double[] objectives;
		double[] variables;
		
		// variables
		variables = new double[solution.numberOfVariables()];
		for (int i = 0; i < solution.numberOfVariables(); i++)
			variables[i] = solution.getDecisionVariables()[i].getValue();
		// objectives
		if (solution.numberOfFeatures() > 0) {
			objectives = new double[solution.numberOfFeatures()];
			for (int i = 0; i < objectives.length; i++)
				objectives[i] = solution.getFeature(i);
		}
		else {
			objectives = new double[solution.numberOfObjectives()];
			for (int i = 0; i < objectives.length; i++)
				objectives[i] = solution.getObjective(i);
		}
		// check if already exists
		if (data_.contains(toHash(variables, objectives)))
			return;
		// append
		for (int i = 0; i < numberOfOutputs_; i++) {
		  inputs_[i].add(variables);
		  outputs_[i].add(objectives[i]);
		  outputsAll_[i].add(objectives[i]);
		}
		inputsAll_.add(variables);
		data_.add(toHash(variables, objectives));
	}
	
	private String toHash(double[] variables, double[] objectives) {
		BASE64Encoder base64;
		StringBuilder concat;
		MessageDigest md;
		
		concat = new StringBuilder();
		base64 = new BASE64Encoder();
		try {
			md = MessageDigest.getInstance("SHA-1");
		}
		catch (NoSuchAlgorithmException e) {
			// tale more obstajat
			md = null;
		}
		for (int i = 0; i < variables.length; i++)
			concat.append(Long.toHexString(Double.doubleToLongBits(variables[i])) + "|");
		concat.append("-|");
		for (int i = 0; i < objectives.length; i++)
			concat.append(Long.toHexString(Double.doubleToLongBits(objectives[i])) + "|");
		return base64.encode(md.digest(concat.toString().getBytes()));
	}
	
	private void optimize(boolean reset) throws Exception {
		Object[] output;
		double[] theta;
		int c = 0;

		//
		// BATCH
		//
		// windowing
  	if (reset) {
  		int start = inputsAll_.size() - windowSize_;
    	start = start >= 0 ? start : 0;
			for (int k = 0; k < numberOfOutputs_; k++) {
			  inputs_[k].clear();
			  outputs_[k].clear();
			  // add last n data (windowSize / numberOfPseudoInputs)
		  	for (int i = start; i < inputsAll_.size(); i++) {
		  		inputs_[k].add(inputsAll_.get(i));
    			outputs_[k].add(outputsAll_[k].get(i));
		  	}
			}
  	}
  	else {
			double[] lfTemp;
			for (int k = 0; k < numberOfOutputs_; k++) {
				c = pseudoOutputs_[k].size();
			  // inputs
				for (int i = 0; i < c; i++) {
					lfTemp = new double[numberOfInputs_];
					for (int j = 0; j < numberOfInputs_; j++) {
						lfTemp[j] = bvs_[k].getDouble(j * c + i + 1);
					}
					inputs_[k].add(lfTemp);
				}
				// outputs
				outputs_[k].addAll(pseudoOutputs_[k]);
				//for (int i = 0; i < pseudoOutputs_.get(0).length; i++)
				//  outputs_[k].add(pseudoOutputs_.get(k)[i]);
		  }
  	}
		/*if (windowSize_ > 0) {
			start = data_.size() - windowSize_;
			start = start >= 0 ? start : 0;
			data_ = data_.subList(start, data_.size());
		}
		else
			data_.clear();*/
		//
		// END OF BATCH
		//
		for (int i = 0; i < numberOfOutputs_; i++) {
			pseudoOutputs_[i] = new ArrayList<Double>();
			c = 0;
			while (true) {
				try {
					c++;
					if (c > 1) {
						output = spgp_.initHyp(1, toArray(outputs_[i]), toArray(inputs_[i]), Math.min(numberOfPseudoInputs_, inputs_[i].size()), Math.min(c * 0.2, 3));
					}
					else {
						if (fixed_ || !resetHyp_)
							output = spgp_.initHyp(1, toArray(outputs_[i]), toArray(inputs_[i]), Math.min(numberOfPseudoInputs_, inputs_[i].size()), hyp_[i]);
						else
							output = spgp_.initHyp(1, toArray(outputs_[i]), toArray(inputs_[i]), Math.min(numberOfPseudoInputs_, inputs_[i].size()));
					}
					theta = ((MWNumericArray)output[0]).getDoubleData();
					theta = minimize(i, theta);
					output = spgp_.unwrap(2, transpose(theta), numberOfInputs_);
					bvs_[i] = (MWNumericArray)output[0];
					if (!fixed_)
						hyp_[i] = (MWNumericArray)output[1];
					break;
				}
				catch (Exception e) {
					//e.printStackTrace();
					e.toString();
				}
			}
  		// GET BV OUTPUTS
			double[] temp;
			//c = bvs_[0].numberOfElements() / numberOfInputs_;
		  //output = spgp_.pred(4, toArray(outputs_[i]), toArray(inputs_[i]), bvs_[i], bvs_[i], hyp_[i], 0, 0);
			output = spgp_.bv(2, toArray(outputs_[i]), toArray(inputs_[i]), bvs_[i], hyp_[i], 0);
			if (((MWCharArray)output[1]).toString().equalsIgnoreCase("OK")) {
				temp = ((MWNumericArray)output[0]).getDoubleData();
				for (int j = 0; j < temp.length; j++)
					pseudoOutputs_[i].add(temp[j]);
			}
			else
       	throw new Exception("to pa ni ok");
		}
	}
		
	private double[] minimize(int index, double[] theta) throws Exception {
		Object[] output;
		double[] w;
		double[] f;

		w = theta.clone();
		output = minimize(w, index);
		w = (double[])output[0];
		f = (double[])output[1];
		System.out.printf("SPGP-%2d: %4.6f --> %4.6f\n", index, f[0], f[f.length-1]);
		return w;
	}
	
	private Object[] minimize(double[] X, int index) throws Exception {
		Object[] output;
		List<Double> fX;
		double INT;
		double EXT;
		int MAX;
		int RATIO;
		double SIG;
		double RHO;
		boolean ls_failed;
		boolean success;
		double[] s;
		double[] X0;
		double[] df0;
		double[] dF0;
		double[] df3;
		double x1;
		double x2;
		double x3;
		double x4;
		double f0;
		double f1;
		double f2;
		double f3;
		double f4;
		double F0;
		double d0;
		double d1;
		double d2;
		double d3;
		double d4;
		double A;
		double B;
		double T;
		int M;
		int i;
	
		d4 = Double.NaN;
		f4 = Double.NaN;
		x4 = Double.NaN;
		fX = new ArrayList<Double>();
		
		INT = 0.1;
		EXT = 3.0;
		MAX = 20;
		RATIO = 10;
		SIG = 0.1;
		RHO = SIG / 2;
		
		i = 0;
		ls_failed = false;
		if (fixed_)
			output = spgp_.lik_nohyp(6, transpose(X), toArray(outputs_[index]), toArray(inputs_[index]), 0, 0);
		else
			output = spgp_.lik(6, transpose(X), toArray(outputs_[index]), toArray(inputs_[index]), 0, 0);
		if (!((MWCharArray)output[5]).toString().equalsIgnoreCase("OK"))
			throw new Exception(((MWCharArray)output[5]).toString());
			
		f0 = ((MWNumericArray)output[0]).getDouble(1);
		//df0 = ((MWNumericArray)output[1]).getDoubleData();
		df0 = concatResults(output);
		fX.add(f0);
		i += (maxOptIter_ < 0 ? 1 : 0);
		s = multiplyVector(df0, -1);
		d0 = multiplyVectors(df0, s);
		x3 = 1 / (1 - d0);
		
		while (i < Math.abs(maxOptIter_)) {
			i += (maxOptIter_ > 0 ? 1 : 0);
			
			X0 = X;
			F0 = f0;
			dF0 = df0;
			if (maxOptIter_ > 0)
				M = MAX;
			else
				M = Math.min(MAX, -maxOptIter_-i);
			
			while (true) {
				x2 = 0;
				f2 = f0;
				d2 = d0;
				f3 = f0;
				df3 = df0;
				success = false;
				while (!success && M > 0) {
					try {
						M -= 1;
						i += (maxOptIter_ < 0 ? 1 : 0);
						//output = spgp_.lik(2, sumVectors(X, multiplyVector(s, x3)), toArray(outputsNoise_, index), toArray(inputs_), bvs_[index].numberOfElements() / numberOfInputs_);
						if (fixed_)
							output = spgp_.lik_nohyp(6, transpose(sumVectors(X, multiplyVector(s, x3))), toArray(outputs_[index]), toArray(inputs_[index]), 0, 0);
						else
							output = spgp_.lik(6, transpose(sumVectors(X, multiplyVector(s, x3))), toArray(outputs_[index]), toArray(inputs_[index]), 0, 0);
						if (!((MWCharArray)output[5]).toString().equals("OK"))
							throw new Exception();
						if (isAnyNanOrInf(((MWNumericArray)output[0]).getDoubleData()) || isAnyNanOrInf(((MWNumericArray)output[1]).getDoubleData()) || isAnyNanOrInf(((MWNumericArray)output[2]).getDoubleData()) || isAnyNanOrInf(((MWNumericArray)output[3]).getDoubleData()) || isAnyNanOrInf(((MWNumericArray)output[4]).getDoubleData()))
							throw new Exception();
						f3 = ((MWNumericArray)output[0]).getDouble(1);
						//df3 = ((MWNumericArray)output[1]).getDoubleData();
						df3 = concatResults(output);
						success = true;
					}
					catch (Exception e) {
						x3 = (x2 + x3) / 2;
					}
				}
				if (f3 < F0) {
					X0 = sumVectors(X, multiplyVector(s, x3));
					F0 = f3;
					dF0 = df3;
				}
				d3 = multiplyVectors(df3, s);  // df3'*s;
				if (d3 > SIG * d0 || f3 > f0 + x3 * RHO * d0 || M == 0)
					break;
				x1 = x2;
				f1 = f2;
				d1 = d2;
				x2 = x3;
				f2 = f3;
				d2 = d3;
				A = 6 * (f1 - f2) + 3 * (d2 + d1) * (x2 - x1);
				B = 3 * (f2 - f1) - (2 * d1 + d2) * (x2 - x1);
				T = B * B - A * d1 * (x2 - x1);
				x3 = x1 - d1 * (x2 - x1) * (x2 - x1) / (B + Math.sqrt(T));
				if (T < 0 || MWNumericArray.isNaN(x3) || MWNumericArray.isInf(x3) || x3 < 0)
					x3 = x2 * EXT;
				else if (x3 > x2 * EXT)
					x3 = x2 * EXT;
				else if (x3 < x2 + INT * (x2 - x1))
					x3 = x2 + INT * (x2 - x1);
			}
			
			while ((Math.abs(d3) > -SIG*d0 || f3 > f0 + x3 * RHO * d0) && M > 0) {
		    if (d3 > 0 || f3 > f0 + x3 * RHO * d0) {
		    	x4 = x3;
		    	f4 = f3;
		    	d4 = d3;
		    }
		    else {
		    	x2 = x3;
		    	f2 = f3;
		    	d2 = d3;
		    }
		    if (f4 > f0)
		    	x3 = x2 - (0.5 * d2 * (x4-x2) * (x4-x2)) / (f4 - f2 - d2 * (x4 - x2));
		    else {
		      A = 6 * (f2 - f4) / (x4 - x2) + 3 * (d4 + d2);
		      B = 3 * (f4 - f2) - (2 * d2 + d4) * (x4 - x2);
		      x3 = x2 + (Math.sqrt(B * B - A * d2 * (x4-x2) * (x4-x2))-B)/A;
		    }
		    if (MWNumericArray.isNaN(x3) || MWNumericArray.isInf(x3))
		    	x3 = (x2 + x4)/2;
		    x3 = Math.max(Math.min(x3, x4 - INT * (x4 - x2)), x2 + INT * (x4 - x2));
		    //output = spgp_.lik(2, sumVectors(X, multiplyVector(s, x3)), toArray(outputsNoise_, index), toArray(inputs_), bvs_[index].numberOfElements() / numberOfInputs_);
	    	if (fixed_)
	    		output = spgp_.lik_nohyp(6, transpose(sumVectors(X, multiplyVector(s, x3))), toArray(outputs_[index]), toArray(inputs_[index]), 0, 0);
	    	else
	    		output = spgp_.lik(6, transpose(sumVectors(X, multiplyVector(s, x3))), toArray(outputs_[index]), toArray(inputs_[index]), 0, 0);
	  		if (!((MWCharArray)output[5]).toString().equalsIgnoreCase("OK"))
	  			throw new Exception(((MWCharArray)output[5]).toString());

	    	f3 = ((MWNumericArray)output[0]).getDouble(1);
				//df3 = ((MWNumericArray)output[1]).getDoubleData();
				df3 = concatResults(output);
				if (f3 < F0) {
		    	X0 = sumVectors(X, multiplyVector(s, x3));
		    	F0 = f3;
		    	dF0 = df3;
		    }
		    M -= 1;
		    i += (maxOptIter_ < 0 ? 1 : 0);
		    d3 = multiplyVectors(df3, s);
			}
			
		  if (Math.abs(d3) < -SIG * d0 && f3 < f0 + x3 * RHO * d0) {
		  	X = sumVectors(X, multiplyVector(s, x3));
		  	f0 = f3;
		  	fX.add(f0);
		  	//fprintf('%s %6i;  Value %4.6e\r', S, i, f0);
		  	//System.out.printf("%s %6d;  Value %4.6f\r", S, i, f0);
		    s = subVectors(multiplyVector(s, (multiplyVectors(df3, df3) - multiplyVectors(df0, df3)) / multiplyVectors(df0, df0)), df3);
		  	df0 = df3;
		  	d3 = d0;
		  	d0 = multiplyVectors(df0, s);
		  	if (d0 > 0) {
		  		s = multiplyVector(df0, -1);
		  		d0 = multiplyVectors(df0, s);
		  	}
	      x3 = x3 * Math.min(RATIO, d3 / (d0 - Double.MIN_VALUE));
	      ls_failed = false;
		  }
	    else {
	    	X = X0;
	    	f0 = F0;
	    	df0 = dF0;
	      if (ls_failed || i > Math.abs(maxOptIter_))
	      	break;
	  		s = multiplyVector(df0, -1);
	  		d0 = multiplyVectors(df0, s);
	      x3 = 1 / (1 - d0);                     
	      ls_failed = true;
	    }
		}
		System.out.printf("SPGP optimization ended after %d/%d iterations\r\n", Math.abs(i), maxOptIter_);
		return new Object[]{X, toArray2(fX), i};
	}
	
	private double[][] transpose(double[] v) {
		double[][] t;
		
		t = new double[v.length][1];
		for (int i = 0; i < v.length; i++)
			t[i][0] = v[i];
		return t;
	}
	
	private double[] concatResults(Object[] res) {
		double[] r;
		int[] ind;
		int[] l;
		int la;
		int c;
		
		c = 0;
		ind = new int[]{1, 2, 3, 4};
		l = new int[ind.length];
		la = 0;
		for (int i = 0; i < l.length; i++) {
		  l[i] = ((MWNumericArray)res[ind[i]]).numberOfElements();
		  la += l[i];
		}
		r = new double[la];
		for (int i = 0; i < ind.length; i++) {
			for (int j = 1; j <= l[i]; j++) {
			  r[c] = ((MWNumericArray)res[ind[i]]).getDouble(j);
			  c++;
			}
		}
		return r;
	}
	
	private boolean isAnyNanOrInf(double[] v) {
		for (int i = 0; i < v.length; i++) {
			if (MWNumericArray.isNaN(v[i]) || MWNumericArray.isInf(v[i]))
				return true;
		}
		return false;
	}
	
	private double multiplyVectors(double[] v1, double[] v2) throws Exception {
		double val;
		
		if (v1.length != v2.length)
			throw new Exception("Vector dimensions must agree.");
		val = 0;
		for (int i = 0; i < v1.length; i++)
			val += v1[i]*v2[i];
		return val;
	}
	
	private double[] multiplyVector(double[] v, double s) {
		double[] m;
		
		m = new double[v.length];
		for (int i = 0; i < v.length; i++)
			m[i] = v[i]*s;
		return m;
	}

	private double[] sumVectors(double[] v1, double[] v2) throws Exception {
		double[] val;
		
		if (v1.length != v2.length)
			throw new Exception("Vector dimensions must agree.");
		val = new double[v1.length];
		for (int i = 0; i < val.length; i++)
			val[i] = v1[i] + v2[i];
		return val;
	}

	private double[] subVectors(double[] v1, double[] v2) throws Exception {
		double[] val;
		
		if (v1.length != v2.length)
			throw new Exception("Vector dimensions must agree.");
		val = new double[v1.length];
		for (int i = 0; i < val.length; i++)
			val[i] = v1[i] - v2[i];
		return val;
	}

	private double[] toArray2(List<Double> list) {
		double[] array;
		
		array = new double[list.size()];
		for (int i = 0; i < array.length; i++)
			array[i] = list.get(i);
		return array;
	}

	private double[][] toArray(ArrayList list) {
		double[][] array;
		
		array = new double[list.size()][];
		if (list.size() > 0) {
			if (list.get(0) instanceof Double) {
				for (int i = 0; i < array.length; i++) {
					array[i] = new double[]{(Double)list.get(i)};
				}
			}
			else if (list.get(0) instanceof double[]) {
				for (int i = 0; i < array.length; i++) {
					array[i] = (double[])list.get(i);
				}				
			}
			else {
				throw new Error("ni pravega tipa");
			}
		}
		return array;
	}

	/*private double[][] toArray(List<double[]> list, int index) {
		double[][] array;
		
		array = new double[list.size()][];
		for (int i = 0; i < array.length; i++)
			array[i] = new double[]{list.get(i)[index]};
		return array;
	}*/

	/*private double[][] toArray(List<Double> list) {
		double[][] array;
		
		array = new double[list.size()][];
		for (int i = 0; i < array.length; i++)
			array[i] = new double[]{list.get(i)};
		return array;
	}*/

	private double[][] toArrayRand(List<double[]> list, int length) {
		double[][] array;
		Random random;
		
		random = new Random();
		array = new double[length][];
		for (int i = 0; i < array.length; i++) {
			array[i] = list.get(random.nextInt(list.size()));
		}
		return array;
	}

	private double[] toArray(Variable[] variables) throws JMException {
		double[] array;
		
		array = new double[variables.length];
		for (int i = 0; i < array.length; i++)
			array[i] = variables[i].getValue();
		return array;
	}
	
	private String toString(List<double[]> list) {
		return toString(list, 0);
	}
	
	private String toString(List<double[]> list, int start) {
		StringBuilder sb;
		double[] values;
		
		sb = new StringBuilder();
		for (int i = start; i < list.size(); i++) {
			values = list.get(i);
			for (int j = 0; j < values.length; j++) {
				sb.append(values[j]);
				sb.append("\t");
			}
			sb.append("\n");
		}
		return sb.toString();
	}
	
} // SPGP