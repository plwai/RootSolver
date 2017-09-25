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
public class NewtonMethod {
    final private int EQNS_LENGTH = 3;
    final private int VARIATE_NUM = 3;
    final private int NEGATE = -1;
    
    public NewtonMethod() {

    }
    
    public double[] f(double[] guess) {
        double [] outVal = new double[VARIATE_NUM];
        
        for(int i = 0; i < outVal.length; i++) {
            outVal[i] = valF(guess, i);
        }
        
        return outVal;
    }
    
    public double fO(double[] guess) {
        double outVal = 0;
        
        for(int i = 0; i < VARIATE_NUM; i++) {
            outVal += valF(guess, i) * valF(guess, i);
        }
        
        return outVal * 0.5;
    }
    
    public double valF(double[] guess, int eqnNum) {
        double result = 0;
        
        switch(eqnNum) {
            case 0:
                result = Math.cos(guess[0]) + Math.cos(guess[1]) + Math.cos(guess[2]) - (3.0/5.0);
                break;
                
            case 1:
                result = Math.cos(3*guess[0]) + Math.cos(3*guess[1]) + Math.cos(3*guess[2]);
                break;
                
            case 2:
                result = Math.cos(5*guess[0]) + Math.cos(5*guess[1]) + Math.cos(5*guess[2]);
                break;
                
            default:
                assert true;
        }
        
        return result * -1;
    }
    
    public double derivative(double [] guess, int row, int col) {
        int location = row * VARIATE_NUM + col;
        double result = 1;
        
        switch(location) {
            case 1:
                result *= Math.sin(guess[0]);
                break;
                
            case 2:
                result *= Math.sin(guess[1]);
                break;
                
            case 3:
                result *= Math.sin(guess[2]);
                break;
            
            case 4:
                result *= (3 * Math.sin(3 * guess[0]));
                break;
                
            case 5:
                result *= (3 * Math.sin(3 * guess[1]));     
                break;
                
            case 6:
                result *= (3 * Math.sin(3 * guess[2]));
                break;
            
            case 7:
                result *= (5 * Math.sin(5 * guess[0]));
                break;
                
            case 8:
                result *= (5 * Math.sin(5 * guess[1]));
                break;
                
            case 9:
                result *= (5 * Math.sin(5 * guess[2]));
                break;
                
            default:
                assert true;
                
        }
        
        return result;
    }
    
    public Matrix jacobian(double [] guess) {
        Matrix matrixJ = new Matrix(EQNS_LENGTH, VARIATE_NUM);
        
        for(int r = 0; r < EQNS_LENGTH; r++) {
            for(int c = 1; c <= VARIATE_NUM; c++) {
                matrixJ.set(r, c-1, (derivative(guess, r, c)));
                //matrixJ.set(r, c-1, (derivative(guess, r, c) ));
            }
        }
        
        return matrixJ;
    }
    
    public Matrix jacobianT(double [] guess) {
        Matrix matrixJ = new Matrix(EQNS_LENGTH, VARIATE_NUM);
        
        for(int r = 0; r < EQNS_LENGTH; r++) {
            for(int c = 1; c <= VARIATE_NUM; c++) {
                matrixJ.set(r, c-1, (derivative(guess, c-1, r+1)));
                //matrixJ.set(r, c-1, (derivative(guess, r, c) ));
            }
        }
        
        return matrixJ;
    }
    
    public double[] solve(double [] guess) {
        double [] fx;
        double [] dgn;
        double [] ans = new double[VARIATE_NUM];
        double damping = 0.1;
        boolean loop = true;
        boolean isSingular = false;
        int itr = 0;
        int stop = 0;
        
        do {
            System.out.println(guess[0]);
            System.out.println(guess[1]);
            System.out.println(guess[2]);
            System.out.println();
            
            Matrix J = jacobian(guess);
            Matrix Jt = jacobianT(guess);
            fx = f(guess);
            
            
            //J.print(3, 10);
            //Jt.print(3, 10);
            //J.lu().getL().print(3, VARIATE_NUM);
            //J.lu().getU().print(3, VARIATE_NUM);
            //(new Matrix(fx, VARIATE_NUM)).print(3, 10);
            

            //dgn = (Jt.times(J).plus(Matrix.identity(EQNS_LENGTH, VARIATE_NUM).times(damping))).lu()
            //            .solve(Jt.times(-1).times(new Matrix(fx, VARIATE_NUM))).getRowPackedCopy();
            
            if(J.lu().isNonsingular() && !isSingular) {
                dgn = J.lu().solve((new Matrix(fx, VARIATE_NUM)).times(NEGATE)).getRowPackedCopy();
                
            }
            else {
                System.out.println("singular");
                isSingular = true;
                dgn = (Jt.times(J).plus(Matrix.identity(EQNS_LENGTH, VARIATE_NUM).times(damping))).lu()
                        .solve(Jt.times(NEGATE).times(new Matrix(fx, VARIATE_NUM))).getRowPackedCopy();

            }

            for(int i = 0; i < VARIATE_NUM; i++) {
                ans[i] = dgn[i] + guess[i];
            }
            
            for(int i = 0; i < VARIATE_NUM; i++) {
                if(Math.abs(dgn[i]) < 0.000000001) {
                    loop = false;
                }
                else {
                    loop = true;
                    break;
                }
            }

            
            System.out.println(dgn[0]);
            System.out.println(dgn[1]);
            System.out.println(dgn[2]);
            System.out.println();
            
            
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
        (new Matrix(fx, VARIATE_NUM)).print(3, 10);
        
        return ans;
    }
    
    public static void main(String[] args) {
        NewtonMethod m = new NewtonMethod();
        
        double [] guess = {1,0.1,0.1};
        double [] ans;
        
        ans = m.solve(guess);
        
        System.out.println(ans[0]);
        System.out.println(ans[1]);
        System.out.println(ans[2]);
    }
}
