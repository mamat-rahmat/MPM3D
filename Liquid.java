package liquid;

import java.awt.event.MouseEvent;
import java.awt.event.KeyEvent;
import java.util.Iterator;
import java.awt.image.ImageObserver;
import java.util.Arrays;
import java.awt.Graphics;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Image;
import java.util.ArrayList;
import java.awt.event.KeyListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseListener;
import javax.swing.JApplet;

public class Liquid extends JApplet implements Runnable, MouseListener, MouseMotionListener, KeyListener {
    ArrayList<Particle> particles;
    int gsizeX;
    int gsizeY;
    int gsizeZ;
    Node[][][] grid;
    ArrayList<Node> active;
    Thread animationThread;
    Image backBuffer;
    Graphics2D backBufferGraphics;
    Material water;
    boolean pressed;
    boolean pressedprev;
    int mx;
    int my;
    int mxprev;
    int myprev;
    Particle[] sorted;
    Color[] speedCol;

    public Liquid() {
        this.particles = new ArrayList<Particle>();
        this.gsizeX = 40;
        this.gsizeY = 40;
        this.gsizeZ = 40;
        this.grid = new Node[this.gsizeX][this.gsizeY][this.gsizeZ];
        this.active = new ArrayList<Node>();
        this.water = new Material(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f);
        this.speedCol = new Color[256];
    }

    public void init() {
        this.setSize(400, 400);
        this.backBuffer = this.createImage(400, 400);
        this.backBufferGraphics = (Graphics2D) this.backBuffer.getGraphics();
        for (int i = 0; i < this.gsizeX; ++i) {
            for (int j = 0; j < this.gsizeY; ++j) {
                for (int k = 0; k < this.gsizeZ; ++k) {
                    this.grid[i][j][k] = new Node();
                }
            }
        }
        (this.animationThread = new Thread(this)).start();
        this.addMouseListener(this);
        this.addMouseMotionListener(this);
        this.addKeyListener(this);
        for (int i = 0; i < 20; ++i) {
            for (int j = 0; j < 25; ++j) {
                for (int k = 0; k < 20; ++k) {
                    Particle p = new Particle(this.water, i + 4, j + 4, k + 4, 0.0f, 0.0f, 0.0f);
                    this.particles.add(p);
                }
            }
        }
        this.sorted = new Particle[this.particles.size()];
        this.particles.toArray(this.sorted);
        for (int i = 0; i < 256; ++i) {
            float frac = i / 256.0f;
            float r = frac;
            float g = 0.5f + 0.5f * frac;
            float b = 1.0f;
            this.speedCol[i] = new Color(r, g, b);
        }
    }

    public void paint(Graphics g) {
        long startTime = System.currentTimeMillis();
        this.simulate();
        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        double fps = 1000.0 / elapsedTime;
        // System.out.println(elapsedTime);

        this.backBufferGraphics.setBackground(Color.black);
        this.backBufferGraphics.clearRect(0, 0, 400, 400);
        this.backBufferGraphics.setColor(Color.BLUE);

        Arrays.sort(this.sorted);
        for (int i = 0; i < this.sorted.length; ++i) {
            Particle p = this.sorted[i];
            float speed = (float) Math.sqrt(p.u * p.u + p.v * p.v + p.w * p.w);
            int index = (int) (speed * 100.0f);
            if (index < 0) {
                index = 0;
            } else if (index > 255) {
                index = 255;
            }
            this.backBufferGraphics.setColor(this.speedCol[index]);
            float vx1 = p.x - 19.5f;
            float vx2 = p.x - p.u - 19.5f;
            float vy1 = p.y - 19.5f;
            float vy2 = p.y - p.v - 19.5f;
            float mul = 10.0f * (1.0f - (p.z - 10.0f) / 100.0f);
            vx1 = 200.0f + vx1 * mul;
            vx2 = 200.0f + vx2 * mul;
            vy1 = 200.0f + vy1 * mul;
            vy2 = 200.0f + vy2 * mul;
            this.backBufferGraphics.drawLine((int) vx1, (int) vy1, (int) vx2, (int) vy2);
        }
        this.backBufferGraphics.setColor(Color.BLACK);
        g.drawImage(this.backBuffer, 0, 0, this);
    }

    public void run() {
        while (true) {
            this.repaint();
        }
    }

    public void simulate() {
        long startTime, stopTime, elapsedTime, time1, time2, time3, time4, time5, time6, time7;
        boolean drag = false;
        float mdx = 0.0f;
        float mdy = 0.0f;
        if (this.pressed && this.pressedprev) {
            drag = true;
            mdx = 0.1f * (this.mx - this.mxprev);
            mdy = 0.1f * (this.my - this.myprev);
        }
        this.pressedprev = this.pressed;
        this.mxprev = this.mx;
        this.myprev = this.my;

        // ----------BEGIN----------
        startTime = System.currentTimeMillis();
        // Reset active node
        for (Node n : this.active) {
            n.az = 0.0f;
            n.ay = 0.0f;
            n.ax = 0.0f;
            n.w = 0.0f;
            n.v = 0.0f;
            n.u = 0.0f;
            n.gz = 0.0f;
            n.gy = 0.0f;
            n.gx = 0.0f;
            n.d = 0.0f;
            n.m = 0.0f;
            n.active = false;
        }
        this.active.clear();
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        time1 = elapsedTime;
        // ---------- END ----------

        // ----------BEGIN----------
        startTime = System.currentTimeMillis();
        // Transfer data from particle to neighbouring node
        // Particle(x, y, z, mat.m) -> Node(m, d, gx, gy, gz)
        for (Particle p : this.particles) {
            p.cx = (int) (p.x - 0.5f);
            p.cy = (int) (p.y - 0.5f);
            p.cz = (int) (p.z - 0.5f);

            float x = p.cx - p.x;
            p.px[0] = 0.5f * x * x + 1.5f * x + 1.125f;
            p.gx[0] = x + 1.5f;
            ++x;
            p.px[1] = -x * x + 0.75f;
            p.gx[1] = -2.0f * x;
            ++x;
            p.px[2] = 0.5f * x * x - 1.5f * x + 1.125f;
            p.gx[2] = x - 1.5f;

            float y = p.cy - p.y;
            p.py[0] = 0.5f * y * y + 1.5f * y + 1.125f;
            p.gy[0] = y + 1.5f;
            ++y;
            p.py[1] = -y * y + 0.75f;
            p.gy[1] = -2.0f * y;
            ++y;
            p.py[2] = 0.5f * y * y - 1.5f * y + 1.125f;
            p.gy[2] = y - 1.5f;

            float z = p.cz - p.z;
            p.pz[0] = 0.5f * z * z + 1.5f * z + 1.125f;
            p.gz[0] = z + 1.5f;
            ++z;
            p.pz[1] = -z * z + 0.75f;
            p.gz[1] = -2.0f * z;
            ++z;
            p.pz[2] = 0.5f * z * z - 1.5f * z + 1.125f;
            p.gz[2] = z - 1.5f;

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        int cxi = p.cx + i;
                        int cyj = p.cy + j;
                        int czk = p.cz + k;
                        Node n = this.grid[cxi][cyj][czk];
                        if (!n.active) {
                            this.active.add(n);
                            n.active = true;
                        }
                        float phi = p.px[i] * p.py[j] * p.pz[k];
                        n.m += phi * p.mat.m;
                        n.d += phi;
                        float dx = p.gx[i] * p.py[j] * p.pz[k];
                        float dy = p.px[i] * p.gy[j] * p.pz[k];
                        float dz = p.px[i] * p.py[j] * p.gz[k];
                        n.gx += dx;
                        n.gy += dy;
                        n.gz += dz;
                    }
                }
            }
            // System.out.format("%.2f %.2f %.2f %d %d %d\n", p.x, p.y, p.z, p.cx, p.cy,
            // p.cz);
        }
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        time2 = elapsedTime;
        // ---------- END ---------

        // ----------BEGIN----------
        startTime = System.currentTimeMillis();
        // Calculate total force
        for (Particle p : this.particles) {
            int cx = (int) p.x;
            int cy = (int) p.y;
            int cz = (int) p.z;
            int cxi = cx + 1;
            int cyi = cy + 1;
            int czi = cz + 1;
            float u = p.x - cx;
            float v = p.y - cy;
            float w = p.z - cz;

            float p000 = this.grid[cx][cy][cz].d;
            float x000 = this.grid[cx][cy][cz].gx;
            float y000 = this.grid[cx][cy][cz].gy;
            float z000 = this.grid[cx][cy][cz].gz;
            float p001 = this.grid[cx][cy][czi].d;
            float x001 = this.grid[cx][cy][czi].gx;
            float y001 = this.grid[cx][cy][czi].gy;
            float z001 = this.grid[cx][cy][czi].gz;

            float p010 = this.grid[cx][cyi][cz].d;
            float x010 = this.grid[cx][cyi][cz].gx;
            float y010 = this.grid[cx][cyi][cz].gy;
            float z010 = this.grid[cx][cyi][cz].gz;
            float p011 = this.grid[cx][cyi][czi].d;
            float x011 = this.grid[cx][cyi][czi].gx;
            float y011 = this.grid[cx][cyi][czi].gy;
            float z011 = this.grid[cx][cyi][czi].gz;

            float p100 = this.grid[cxi][cy][cz].d;
            float x100 = this.grid[cxi][cy][cz].gx;
            float y100 = this.grid[cxi][cy][cz].gy;
            float z100 = this.grid[cxi][cy][cz].gz;
            float p101 = this.grid[cxi][cy][czi].d;
            float x101 = this.grid[cxi][cy][czi].gx;
            float y101 = this.grid[cxi][cy][czi].gy;
            float z101 = this.grid[cxi][cy][czi].gz;

            float p110 = this.grid[cxi][cyi][cz].d;
            float x110 = this.grid[cxi][cyi][cz].gx;
            float y110 = this.grid[cxi][cyi][cz].gy;
            float z110 = this.grid[cxi][cyi][cz].gz;
            float p111 = this.grid[cxi][cyi][czi].d;
            float x111 = this.grid[cxi][cyi][czi].gx;
            float y111 = this.grid[cxi][cyi][czi].gy;
            float z111 = this.grid[cxi][cyi][czi].gz;

            float dx00 = p100 - p000;
            float dx01 = p101 - p001;
            float dx10 = p110 - p010;

            float dy00 = p010 - p000;
            float dy01 = p011 - p001;
            float dy10 = p110 - p100;

            float dz00 = p001 - p000;
            float dz01 = p011 - p010;
            float dz10 = p101 - p100;

            float C000 = p000;
            float C100 = x000;
            float C010 = y000;
            float C001 = z000;

            float C310 = x110 - x100 + x010 - x000 - 2.0F * (dx10 - dx00);
            float C210 = 3.0F * (dx10 - dx00) - 2.0F * (x010 - x000) - (x110 - x100);

            float C301 = x101 - x100 + x001 - x000 - 2.0F * (dx01 - dx00);
            float C201 = 3.0F * (dx01 - dx00) - 2.0F * (x001 - x000) - (x101 - x100);

            float C130 = y110 - y010 + y100 - y000 - 2.0F * (dy10 - dy00);
            float C120 = 3.0F * (dy10 - dy00) - 2.0F * (y100 - y000) - (y110 - y010);

            float C031 = y011 - y010 + y001 - y000 - 2.0F * (dy01 - dy00);
            float C021 = 3.0F * (dy01 - dy00) - 2.0F * (y001 - y000) - (y011 - y010);

            float C103 = z101 - z001 + z100 - z000 - 2.0F * (dz10 - dz00);
            float C102 = 3.0F * (dz10 - dz00) - 2.0F * (z100 - z000) - (z101 - z001);

            float C013 = z011 - z001 + z010 - z000 - 2.0F * (dz01 - dz00);
            float C012 = 3.0F * (dz01 - dz00) - 2.0F * (z010 - z000) - (z011 - z001);

            float C300 = x100 + x000 - 2.0F * dx00;
            float C200 = 3.0F * dx00 - x100 - 2.0F * x000;

            float C030 = y010 + y000 - 2.0F * dy00;
            float C020 = 3.0F * dy00 - y010 - 2.0F * y000;

            float C003 = z001 + z000 - 2.0F * dz00;
            float C002 = 3.0F * dz00 - z001 - 2.0F * z000;

            float C110 = x010 - C100 - C120 - C130;
            float C011 = y001 - C010 - C012 - C013;
            float C101 = z100 - C001 - C201 - C301;

            float A = p100 + y100 + z100 + C011 + C020 + C002 + C120 + C021 + C102 + C012 + C030 + C003 + C130 + C031
                    + C103 + C013;

            float f111_A = p111 - A;

            float x0 = x111 - x110 - x101 + x100;
            float x1 = x011 - x010 - x001 + x000;
            float C311 = x0 + x1 - 2.0F * f111_A;
            float C211 = 3.0F * f111_A - x0 - 2.0F * x1;

            float y0 = y111 - y110 - y011 + y010;
            float y1 = y101 - y100 - y001 + y000;
            float C131 = y0 + y1 - 2.0F * f111_A;
            float C121 = 3.0F * f111_A - y0 - 2.0F * y1;

            float z0 = z111 - z101 - z011 + z001;
            float z1 = z110 - z100 - z010 + z000;
            float C113 = z0 + z1 - 2.0F * f111_A;
            float C112 = 3.0F * f111_A - z0 - 2.0F * z1;

            float C111 = x1 + y1 + z1 - 2.0F * f111_A;

            float density = C000 + (C001 + (C002 + C003 * w) * w) * w + (C010 + (C011 + (C012 + C013 * w) * w) * w) * v
                    + (C020 + C021 * w) * v * v + (C030 + C031 * w) * v * v * v
                    + (C100 + (C110 + (C120 + C130 * v) * v) * v + (C101 + (C111 + (C121 + C131 * v) * v) * v) * w
                            + (C102 + C112 * v) * w * w + (C103 + C113 * v) * w * w * w) * u
                    + (C200 + C210 * v + (C201 + C211 * v) * w) * u * u
                    + (C300 + C310 * v + (C301 + C311 * v) * w) * u * u * u;

            float pressure = p.mat.m / p.mat.rd * p.mat.k * (density - p.mat.rd);
            float fx = 0.0f;
            float fy = 0.0f;
            float fz = 0.0f;

            // Rebound or reflect after striking a surface
            if (p.x < 3.0f) {
                fx += p.mat.m * (3.0f - p.x);
            } else if (p.x > this.gsizeX - 4) {
                fx += p.mat.m * (this.gsizeX - 4 - p.x);
            }
            if (p.y < 3.0f) {
                fy += p.mat.m * (3.0f - p.y);
            } else if (p.y > this.gsizeY - 4) {
                fy += p.mat.m * (this.gsizeY - 4 - p.y);
            }
            if (p.z < 3.0f) {
                fz += p.mat.m * (3.0f - p.z);
            } else if (p.z > this.gsizeZ - 4) {
                fz += p.mat.m * (this.gsizeZ - 4 - p.z);
            }

            // Calculate external force
            if (drag) {
                float vx = Math.abs(p.x - 0.1f * this.mx);
                float vy = Math.abs(p.y - 0.1f * this.my);
                float vz = Math.abs(p.z - 20.0f);
                if (vx < 10.0f && vy < 10.0f && vz < 10.0f) {
                    float weight = p.mat.m * (1.0f - vx / 10.0f) * (1.0f - vy / 10.0f) * (1.0f - vz / 10.0f);
                    fx += weight * (mdx - p.u);
                    fy += weight * (mdy - p.v);
                }
            }

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        Node n = this.grid[p.cx + i][p.cy + j][p.cz + k];
                        float phi = p.px[i] * p.py[j] * p.pz[k];
                        float gx = p.gx[i] * p.py[j] * p.pz[k];
                        float gy = p.px[i] * p.gy[j] * p.pz[k];
                        float gz = p.px[i] * p.py[j] * p.gz[k];
                        // Total force
                        n.ax += -(gx * pressure) + fx * phi;
                        n.ay += -(gy * pressure) + fy * phi;
                        n.az += -(gz * pressure) + fz * phi;
                    }
                }
            }
        }
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        time3 = elapsedTime;
        // ---------- END ---------

        // ----------BEGIN----------
        startTime = System.currentTimeMillis();
        // Calculate node accelaration = (total force) / mass
        Iterator<Node> it = this.active.iterator();
        while (it.hasNext()) {
            Node n = it.next();
            if (n.m > 0.0f) {
                n.ax /= n.m;
                n.ay /= n.m;
                n.az /= n.m;
                n.ay += 0.1f; // gravity
            }
        }
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        time4 = elapsedTime;
        // ---------- END ---------

        // ----------BEGIN----------
        startTime = System.currentTimeMillis();
        // Calculate particle velocity, particle momentum and node ???
        for (Particle p : this.particles) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        Node n = this.grid[p.cx + i][p.cy + j][p.cz + k];
                        float phi = p.px[i] * p.py[j] * p.pz[k];
                        // Particle velocity
                        p.u += phi * n.ax;
                        p.v += phi * n.ay;
                        p.w += phi * n.az;
                    }
                }
            }
            // Particle momentum
            float mu = p.mat.m * p.u;
            float mv = p.mat.m * p.v;
            float mw = p.mat.m * p.w;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        Node n = this.grid[p.cx + i][p.cy + j][p.cz + k];
                        float phi = p.px[i] * p.py[j] * p.pz[k];
                        // Node velocity
                        n.u += phi * mu;
                        n.v += phi * mv;
                        n.w += phi * mw;
                    }
                }
            }
        }
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        time5 = elapsedTime;
        // ---------- END ---------

        // ----------BEGIN----------
        startTime = System.currentTimeMillis();
        // Calculate node velocity
        it = this.active.iterator();
        while (it.hasNext()) {
            Node n = it.next();
            if (n.m > 0.0f) {
                n.u /= n.m;
                n.v /= n.m;
                n.w /= n.m;
            }
        }
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        time6 = elapsedTime;
        // ---------- END ---------

        // ----------BEGIN----------
        startTime = System.currentTimeMillis();
        // Calculate particle position
        // Node velocity -> Particle velocity -> Particle Position
        for (Particle p : this.particles) {
            float gu = 0.0f;
            float gv = 0.0f;
            float gw = 0.0f;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        Node n = this.grid[p.cx + i][p.cy + j][p.cz + k];
                        float phi = p.px[i] * p.py[j] * p.pz[k];
                        // Particle velocity
                        gu += phi * n.u;
                        gv += phi * n.v;
                        gw += phi * n.w;
                    }
                }
            }
            // Particle Position
            p.x += gu;
            p.y += gv;
            p.z += gw;
            p.u += 1.0f * (gu - p.u);
            p.v += 1.0f * (gv - p.v);
            p.w += 1.0f * (gw - p.w);

            // If particle move outsite box, set to border
            if (p.x < 1.0f) {
                p.x = 1.0f + (float) Math.random() * 0.01f;
                p.u = 0.0f;
            } else if (p.x > this.gsizeX - 2) {
                p.x = this.gsizeX - 2 - (float) Math.random() * 0.01f;
                p.u = 0.0f;
            }
            if (p.y < 1.0f) {
                p.y = 1.0f + (float) Math.random() * 0.01f;
                p.v = 0.0f;
            } else if (p.y > this.gsizeY - 2) {
                p.y = this.gsizeY - 2 - (float) Math.random() * 0.01f;
                p.v = 0.0f;
            }
            if (p.z < 1.0f) {
                p.z = 1.0f + (float) Math.random() * 0.01f;
                p.w = 0.0f;
            } else if (p.z > this.gsizeZ - 2) {
                p.z = this.gsizeZ - 2 - (float) Math.random() * 0.01f;
                p.w = 0.0f;
            }
        }
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        time7 = elapsedTime;
        // ---------- END ---------

        // System.out.format("%d\t%d\t%d\t%d\t%d\t%d\t%d\t", time1, time2, time3, time4,
        // time5, time6, time7);
    }

    public void keyPressed(KeyEvent arg0) {
    }

    public void keyReleased(KeyEvent arg0) {
    }

    public void keyTyped(KeyEvent arg0) {
    }

    public void mouseDragged(MouseEvent arg0) {
        this.pressed = true;
        this.mx = arg0.getX();
        this.my = arg0.getY();
    }

    public void mouseMoved(MouseEvent arg0) {
        this.mx = arg0.getX();
        this.my = arg0.getY();
    }

    public void mouseClicked(MouseEvent arg0) {
    }

    public void mouseEntered(MouseEvent arg0) {
    }

    public void mouseExited(MouseEvent arg0) {
    }

    public void mousePressed(MouseEvent arg0) {
        this.pressed = true;
    }

    public void mouseReleased(MouseEvent arg0) {
        this.pressed = false;
    }
}
