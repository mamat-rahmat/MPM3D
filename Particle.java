package liquid;

import java.awt.Color;


public class Particle implements Comparable<Particle>
{
    public Material mat;
    public float x;
    public float y;
    public float z;
    public float u;
    public float v;
    public float w;
    public int cx;
    public int cy;
    public int cz;
    public float[] px;
    public float[] py;
    public float[] pz;
    public float[] gx;
    public float[] gy;
    public float[] gz;
    public Color c;
    
    public Particle(final Material mat, final float x, final float y, final float z, final float u, final float v, final float w) {
        this.px = new float[3];
        this.py = new float[3];
        this.pz = new float[3];
        this.gx = new float[3];
        this.gy = new float[3];
        this.gz = new float[3];
        this.mat = mat;
        this.x = x;
        this.y = y;
        this.z = z;
        this.u = u;
        this.v = v;
        this.w = w;
        this.c = Color.getHSBColor((z - 4.0f) / 10.0f, 1.0f, 1.0f);
    }
    
    public int compareTo(final Particle p2) {
        return (this.z < p2.z) ? 1 : 0;
    }
}
