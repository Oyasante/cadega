/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package padll;

import java.util.ArrayList;

/**
 *
 * @author Yaw's PC
 */
class GeneLink {
    private String geneName;
    private ArrayList<TissueLink> tissues;
    private ArrayList<String> tList;
    private ArrayList<Double> vals;
    private double sum;
    private int count;
    private double average;
            
    public GeneLink(String name) {
        geneName = name;
        vals = new ArrayList<>();
        tList = new ArrayList<>();
    }
    
    public String getGene() {
        return geneName;
    }
    
    public void addTissue(String t, double v) {
        tissues.add(new TissueLink(t, v));
    }
    
    public void setTissue(String t) {
        tList.add(t);
    }
    
    public void setVal(double d) {
        vals.add(d);
        sum += d;
        count++;
        average = sum / count;
        
    }
    
    public double getSum() {
        return sum;
    }
    
    public double getAvg() {
        return average;
    }
    
    public double getInVal(int a) {
        return vals.get(a);
    }
    
    public boolean hasTissue(String t) {
        for (String tis: tList) {
            if (tis.equals(t)) {
                return true;
            }
        }
        return false;
    }
    
    public boolean hasTissue(TissueLink t) {
        boolean ret = false;
        for (TissueLink current: tissues) {
            String str = current.getTissue();
            if ((t.getTissue()).equals(str)) {
                ret = true;
                break;
            }
        }
        return ret;
    }
    
    public double getSumValue() {
        double ret = 0;
        for (TissueLink current: tissues) {
            ret += current.getSensitivity();
        }
        return ret;
    }
    public double getAverageValue() {
        double ret = 0;
        ret = this.getSumValue() / tissues.size();
        
        return ret;
    }
}
