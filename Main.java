/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package padll;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

/**
 *
 * @author Yaw's PC
 */
public class Main {
    
    public static final String SIGIL = "1";
    /*     */
    public static final String PATH = "C:\\Users\\Yaw's PC\\Documents\\NETPHIX\\data\\";
    public static final String PATH2 = "C:\\Users\\Yaw's PC\\Documents\\NETPHIX\\";
    public static final String PATH4 = "C:\\Users\\Yaw's PC\\Documents\\CADEGA\\inputs\\";
    public static final String PATH5 = "C:\\Users\\Yaw's PC\\Documents\\CADEGA\\drugs\\";
    
    public static final String ALTERATION = "AlterationsV2_final.txt\\AlterationsV2_clean.txt";
    public static final String[] DRUGS = {"Drug01","Drug02","Drug03","Drug04","Drug05","Drug06"};
    // public static final String[] DRUGS2 = {"Drug01","Drug02","Drug03","Drug04","Drug05","Drug06","Drug07","Drug08","Drug09","Drug10"};
    public static final String[] DRUGNAME = {"(5Z)-7-Oxozeaenol", "5-Fluorouracil", "681640", "A-443654","A-770041","Afatinib"};
    
    private static ArrayList<String> drugTarget;

    private static String padllTrans(ArrayList<GeneLink> genes) {
        StringBuilder ret = new StringBuilder();
        
        // set up inverse scaled gaussian of the scaled values
        // (1/(stdev*sqrt(2pi))*e^(-.5*Z)
        double max = Double.MIN_VALUE;
        double min = Double.MAX_VALUE;
        
        double std = 0;
        double sum = 0;
        int nonZeroCount = 0;
        ArrayList<Double> gauss = new ArrayList<>();
        
        ArrayList<Double> scale = new ArrayList<>();
        double premax = Double.MIN_VALUE;
        
        for (GeneLink gh: genes) {
            if (premax < gh.getSum())
                premax = gh.getSum();
            scale.add(gh.getSum());
        }
        
        for (Double gh: scale) {
            double val = gh / premax;
            if (Double.isNaN(val))
                continue;
            sum += val;
            if (val != 0)
                nonZeroCount++;
        }
        
        double mean = sum/nonZeroCount;
        double stdev = 0;
        ArrayList<String> geneGet = new ArrayList<>();
        int count = 0;
        
        for (Double gh: scale) {
            double val = gh / premax;
            if (Double.isNaN(val)) {
                count++;
                continue;
            }
            
            double X;
            if (val != 0) {
                X = (val - mean) * (val - mean);
                geneGet.add(genes.get(count).getGene());
                gauss.add(X);
                stdev += X; 
            }
            count++;
        }
        
        stdev = Math.sqrt(stdev / (nonZeroCount));
        double term1 = (1/(stdev*Math.sqrt(2*Math.PI)));
        ArrayList<Double> output = new ArrayList<>();
        
        for (int j = 0; j < gauss.size(); j++) {
            double d = (term1 * Math.exp(-.5*gauss.get(j)/(stdev*stdev)));
            if (min > d)
                min = d;
            if (max < d)
                max = d;
            output.add(d);
        }
        
        for (int k = 0; k < output.size();k++) {
            ret.append(geneGet.get(k) + "\t" + (1-((output.get(k) - min)/(max-min))) + "\n");
        }
        
        
        
        return ret.toString();
    }
    private ArrayList<Gene> all;
    private static ArrayList<GeneLink> genes;
    private static ArrayList<TissueLink> tissueToDrug;
    
    public static void main(String[] args) throws FileNotFoundException, IOException {
        // write the final drug data set
        
        tissueToDrug = new ArrayList<>();
        genes = new ArrayList<>();
        drugTarget = new ArrayList<>();
        
        for(int k = 0;k < DRUGS.length; k++) {
            drugTarget.add("CDTarget" + k);
        }
        System.out.println("Printed " + drugTarget.size() + " drugs.");
        
        // test condition
        boolean onePass = false;
        int count = 0;
        scanAlteration();
        System.out.println("Scanned alterations");
        
        for (String target: drugTarget) {
            //if (count != 5) continue;
            scanDrugToTissues(target);
            System.out.println("Ran drug to tissues");
            

            for (TissueLink first: tissueToDrug) {
                String tis = first.getTissue();
                double val = first.getSensitivity();

                for (GeneLink g: genes) {
                    if(g.hasTissue(tis)) {
                        g.setVal(val);
                    }
                }
            }
            
            String rows = padllTrans(genes);
            //String rows2 = paddlTrans2(genes);

                
            // for netphix File out = new File(PATH2+ "output" + count + ".txt");
            File out = new File(PATH4+ "CDInput" + count + ".tsv");
            FileWriter fw = new FileWriter(out);
            fw.write(rows);
            /*
            //fw.write("GENE\tSCORE\n");
            for (GeneLink g: genes) {
                // for testing purposes fw.write(g.getGene() + "\t" + g.getSum() + "\t" + g.getAvg() +"\n");
                
                // scores are the normalized sums of AUC values
                //double outputVal = ((g.getAvg() - total/newScores.length)/std);
                //outputVal = outputVal * outputVal; //(outputVal + Math.abs(minVal))/ (maxVal + + Math.abs(minVal));
                
                double outputVal2 = g.getSum() / max;
                
                // writes to the output file
                if (!Double.isNaN(outputVal2) && outputVal2 != 0) //remember to add anti-zeros back in*****
                    fw.write(g.getGene() + "\t" + (outputVal2) + "\n");
                    //fw.write(g.getGene() + "\t" + (1/(std*Math.sqrt(2*Math.PI)))*Math.exp(-.5*outputVal) + "\n");
                fw.flush();
            }
            */
            fw.flush();
            fw.close();
            
            count++;
            if (onePass) break;
        }
    }
    
    public static void scanDrugToTissues(String fileName) throws FileNotFoundException {
        String path = "gdsc_auc\\";
        // FOR netphix File f = new File(PATH+ path + fileName + ".txt");
        File f = new File(PATH5 + fileName + ".txt");
        Scanner sc = new Scanner(f);
        sc.useDelimiter("\n");
        int outerLoopCount = 0;
        
        while (sc.hasNext()) {
            Scanner sc2 = new Scanner(sc.next());
            sc2.useDelimiter("\t");
            int innerLoopCount = 0;
            
            while (sc2.hasNext()) {
                String current = sc2.next();
                
                if (innerLoopCount == 0 && outerLoopCount == 0) {
                    innerLoopCount++;
                    continue;
                } else if (innerLoopCount == 0 && outerLoopCount > 0){
                    for (int i = 0; i < tissueToDrug.size();i++) {
                        tissueToDrug.get(i).setDrug(current);
                    }
                } else if (innerLoopCount > 0 && outerLoopCount == 0) {
                    tissueToDrug.add(new TissueLink(current));
                } else if (innerLoopCount > 0 && outerLoopCount > 0) {
                    tissueToDrug.get(innerLoopCount-1).setValue(Double.parseDouble(current));
                }
                innerLoopCount++;
            }
            outerLoopCount++;
        }
    }
    
    public static void scanAlteration () throws FileNotFoundException {
        File f = new File(PATH+ ALTERATION);
        Scanner sc = new Scanner(f);
        sc.useDelimiter("\n");
        int outerLoopCount = 0;
        ArrayList<String> tis = new ArrayList<>();
        
        while (sc.hasNext()) {
            Scanner sc2 = new Scanner(sc.next());
            sc2.useDelimiter("\t");
            int innerLoopCount = 0;
            
            while (sc2.hasNext()) {
                String current = sc2.next();
                
                if (innerLoopCount == 0 && outerLoopCount == 0) {
                    innerLoopCount++;
                    continue;
                } else if (innerLoopCount == 0 && outerLoopCount > 0){
                    genes.add(new GeneLink(current));
                } else if (innerLoopCount > 0 && outerLoopCount == 0) {
                    tis.add(current);
                } else if (innerLoopCount > 0 && outerLoopCount > 0) {
                    int value = (int)Double.parseDouble(current);
                    
                    if (value > 0) {
                        genes.get(outerLoopCount-1).setTissue(tis.get(innerLoopCount-1));
                    }
                    
                }
                innerLoopCount++;
            }
            outerLoopCount++;
        }
    }
    
    
}
