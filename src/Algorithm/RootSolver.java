/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Algorithm;

import Jama.Matrix;

/**
 *
 * @author Wai Pai Lee
 */
public class RootSolver {
    final private int NEGATE = -1;
    
    public RootSolver() {}
    
    public double calMagnitudeSq(double [] fx) {
        double result = 0;
        
        for(int i = 0; i < fx.length; i++) {
            result += (fx[i] * fx[i]);
        }
        
        return result;
    }
    
    /* Solve provided equation that implemented Equation interface
     * @param eqn, Equation class with all required method implemented
     * @param guess, starting guess value. 0 is not allowed.
     * @param algo, the chosen algorithm. Currently only newton method and LM method are available.
     *
     * @return solution of the equations
     */
    public double[] solve(Equations eqn, double [] guess, int algo) {
        double [] fx;
        double [] dgn = new double[eqn.getVariateNum()];
        double [] ans = new double[eqn.getVariateNum()];
        double damping = 0.01;
        double [] squareError;
        double [] error;
        double preError = 0;
        boolean loop = true;
        boolean isSingular = false;
        int itr = 0;
        int stop = 0;
        
        do {
            Matrix J = eqn.jacobian(guess);
            Matrix Jt = eqn.jacobianT(guess);
            fx = eqn.f(guess);
            
            switch(algo) {
                case 1:
                    // Newton Raphson With Levenberg-Marquardt Method
                    if(J.lu().isNonsingular() && !isSingular) {
                        dgn = J.lu().solve((new Matrix(fx, eqn.getVariateNum())).times(NEGATE)).getRowPackedCopy();

                    }
                    else {
                        isSingular = true;
                        dgn = (Jt.times(J).plus(Matrix.identity(eqn.getEqnNum(), eqn.getVariateNum()).times(damping))).lu()
                                .solve(Jt.times(NEGATE).times(new Matrix(fx, eqn.getVariateNum()))).getRowPackedCopy();

                    }
                    
                    break;
                    
                case 2:
                    // Levenberg-Marquardt Method
                    dgn = (Jt.times(J).plus(Matrix.identity(eqn.getEqnNum(), eqn.getVariateNum()).times(damping))).lu()
                        .solve(Jt.times(-1).times(new Matrix(fx, eqn.getVariateNum()))).getRowPackedCopy();
                    
                    break;
                    
                default:
                    // Error: Algorithm is not supported
                    return dgn;
            }

            for(int i = 0; i < eqn.getVariateNum(); i++) {
                ans[i] = dgn[i] + guess[i];
            }
            
            error = new Matrix(eqn.getVariateNum(), eqn.getEqnNum()).getRowPackedCopy();
            
            Matrix Jd = eqn.jacobian(dgn);

            for(int r = 0; r < eqn.getEqnNum(); r++) {
                int [] row = {r};
                
                squareError = ((Jd.getMatrix(row, 0, J.getColumnDimension() - 1).transpose()
                        .times(Jd.getMatrix(row, 0, J.getColumnDimension() - 1)))
                        .minus((new Matrix(fx, eqn.getVariateNum())).times(-2)
                        .times(Jd.getMatrix(row, 0, J.getColumnDimension() - 1))))
                        .plus(Matrix.identity(eqn.getEqnNum(), eqn.getVariateNum()).times(calMagnitudeSq(guess)))
                        .getRowPackedCopy();
                        
                error = (new Matrix(squareError, eqn.getVariateNum()))
                        .plus(new Matrix(error, eqn.getVariateNum()))
                        .getRowPackedCopy();
            }
            
            error = (new Matrix(error, eqn.getVariateNum())).times(1.0/eqn.getVariateNum()).getColumnPackedCopy();
            
            System.out.println(Math.abs(preError - calMagnitudeSq(error)));
            System.out.println();
            
            if(Math.abs(preError - calMagnitudeSq(error)) < 10) {
                damping = 1.0;
            }
            else {
                damping = 0.01;
            }
            
            preError = calMagnitudeSq(error);
            
            for(int i = 0; i < eqn.getVariateNum(); i++) {
                if(Math.abs(dgn[i]) < 0.000000001) {
                    loop = false;
                }
                else {
                    loop = true;
                    break;
                }
            }
            
            guess = ans.clone();
            
            if(!loop) {
                stop++;
                loop = true;
            }
            else {
                stop = 0;
            }
            
            if(stop == 4) {
                break;
            }
            
            itr++;
        }while(loop);
                
        System.out.println(itr);
        (new Matrix(fx, eqn.getVariateNum())).print(3, 10);
        
        return ans;
    }
}
