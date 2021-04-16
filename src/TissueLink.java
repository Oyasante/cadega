/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package padll;

/**
 *
 * @author Yaw's PC
 */
class TissueLink {
    private String tissue;
    private String drug;
    private double dsValue;
    
    public TissueLink() {
        
    }
    
    public TissueLink(String tis) {
        tissue = tis;
    }
    
    public TissueLink(String tis, double val) {
        tissue = tis;
        dsValue = val;
    }
    
    
    public void setValue(double ds) {
        dsValue = ds;
    }
    
    public void setDrug(String dd) {
        drug = dd;
    }
    
    public void setTissue(String tis) {
        tissue = tis;
    }
    
    
    public String getDrug () {
        return drug;
    }
    
    public String getTissue() {
        return tissue;
    }
    
   
    public double getSensitivity() {
        return dsValue;
    }
}
