import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.security.*;
import java.util.*;
import javax.imageio.*;

// --------------------------------------------------------
class Pnt {
    public int x,y;
    public Pnt() {};
    public Pnt(int x1, int y1) {
        x = x1;
        y = y1;
    }
    public boolean equals(Pnt other) {
        return (x == other.x && y == other.y);
    }
    public int dist2(Pnt other) {
        return (x-other.x)*(x-other.x) + (y-other.y)*(y-other.y);
    }
    public String toString() {
        return "(" + x + "," + y + ")";
    }
}

public class RoadsAndJunctionsVis {
    static int minS = 100, maxS = 1000;
    static int minNC = 10, maxNC = 100;
    static double minJC = 0, maxJC = 10;
    static double minFailProb = 0, maxFailProb = 0.4;

    int S;                  // size of the area
    double jCost;           // cost of building a junction
    double jFailProb;       // probability of a junction construction failing
    int NC;                 // number of cities
    Pnt[] cities;           // coordinates of cities
    int NJ;                 // number of junctions
    Pnt[] junct;            // coordinates of newly built junctions
    HashSet<Integer> cUsed; // coordinates of cities and junctions as a set
    int[] junctOk;          // indicates whether the junction has been built successfully and can be used
    Pnt[] vertices;         // cities and junctions stored together
    SecureRandom rnd;

    int NR;                 // number of roads built
    int[] roadS, roadE;     // start and end cities/junctions of newly built roads
    // -----------------------------------------
    boolean isInside(int x, int y) {
        return (x >= 0 && x <= S && x >= 0 && x <= S);
    }
    int cityToInt(Pnt p) {
        return p.x * (S+1) + p.y;
    }
    int roadToInt(int start, int end) {
        return start * (NC + NJ) + end;
    }
    // -----------------------------------------
    String generate(String seedStr) {
      try {
        // generate test case
        rnd = SecureRandom.getInstance("SHA1PRNG");
        long seed = Long.parseLong(seedStr);
        rnd.setSeed(seed);
        S = rnd.nextInt(maxS - minS + 1) + minS;
        NC = rnd.nextInt(maxNC - minNC + 1) + minNC;
        jCost = rnd.nextDouble() * (maxJC - minJC) + minJC;
        jFailProb = rnd.nextDouble() * (maxFailProb - minFailProb) + minFailProb;
        if (seed == 1) {
            S = minS;
            NC = minNC;
            jCost = minJC;
            jFailProb = minFailProb;
        }
        else if (seed == 2) {
            S = maxS;
            NC = maxNC;
            jCost = maxJC;
            jFailProb = maxFailProb;
        }

        // generate the cities' coordinates (all distinct)
        cities = new Pnt[NC];
        cUsed = new HashSet<>();
        for (int i = 0; i < NC; ++i) {
            while (true) {
                cities[i] = new Pnt(rnd.nextInt(S + 1), rnd.nextInt(S + 1));
                Integer point = new Integer(cityToInt(cities[i]));
                if (cUsed.contains(point))
                    continue;
                cUsed.add(point);
                break;
            }
        }

        StringBuffer sb = new StringBuffer();
        sb.append("S = ").append(S).append('\n');
        sb.append("Number of cities = ").append(NC).append('\n');
        sb.append("Junction cost = ").append(jCost).append('\n');
        sb.append("Probability of junction construction failure = ").append(jFailProb).append('\n');
        sb.append("Cities:\n");
        for (int i = 0; i < NC; ++i)
            sb.append(cities[i] + " ");
        sb.append('\n');
        return sb.toString();
      }
      catch (Exception e) {
        addFatalError("An exception occurred while generating test case.");
        e.printStackTrace(); 
        return "";
      }
    }
    // -----------------------------------------
    public double runTest(String seed) {
      try {
        String test = generate(seed);
        if (debug)
            System.out.println(test);

        if (proc != null) {
            // call the solution
            int[] citiesArg = new int[2 * NC];
            for (int i = 0; i < NC; ++i) {
                citiesArg[2*i] = cities[i].x;
                citiesArg[2*i+1] = cities[i].y;
            }

            // first build the junctions
            int[] junctionsRet;
            try { 
                junctionsRet = buildJunctions(S, citiesArg, jCost, jFailProb);
            } catch (Exception e) {
                addFatalError("Failed to get result from buildJunctions.");
                return -1.0;
            }

            // check the return and convert it to junctions
            if (junctionsRet == null) {
                addFatalError("Your return from buildJunctions contained invalid number of elements.");
                return -1.0;
            }
            if (junctionsRet.length % 2 == 1) {
                addFatalError("Your return from buildJunctions contained odd number of elements.");
                return -1.0;
            }
            NJ = junctionsRet.length / 2;
            if (NJ > 2 * NC) {
                addFatalError("You can build at most " + 2 * NC + " new junctions.");
                return -1.0;
            }
            junct = new Pnt[NJ];
            for (int i = 0; i < NJ; ++i) {
                if (!isInside(junctionsRet[2*i], junctionsRet[2*i+1])) {
                    addFatalError("You can only build junctions inside the area.");
                    return -1.0;
                }
                junct[i] = new Pnt(junctionsRet[2*i], junctionsRet[2*i+1]);
                Integer point = new Integer(cityToInt(junct[i]));
                if (cUsed.contains(point)) {
                    addFatalError("You can only build junctions on places not occupied by cities or other junctions.");
                    return -1.0;
                }
                cUsed.add(point);
            }
            // decide which junctions are functional and which are not
            junctOk = new int[NJ];
            for (int i = 0; i < NJ; ++i)
                junctOk[i] = (rnd.nextDouble() >= jFailProb) ? 1 : 0;

            // next build the roads (passing information about junctions to the method)
            int[] roadsRet;
            try { 
                roadsRet = buildRoads(junctOk);
            } catch (Exception e) {
                addFatalError("Failed to get result from buildRoads.");
                return -1.0;
            }
            // check the return and convert it to junctions
            if (roadsRet == null) {
                addFatalError("Your return from buildRoads contained invalid number of elements.");
                return -1.0;
            }
            if (roadsRet.length % 2 == 1) {
                addFatalError("Your return from buildRoads contained odd number of elements.");
                return -1.0;
            }
            NR = roadsRet.length / 2;
            int maxR = (NC + NJ) * (NC + NJ - 1) / 2;
            if (NR > maxR) {
                addFatalError("You can build at most " + maxR + " roads.");
                return -1.0;
            }
            
            HashSet<Integer> rUsed = new HashSet<>();
            roadS = new int[NR];
            roadE = new int[NR];
            for (int i = 0; i < NR; ++i) {
                if (roadsRet[2*i] < 0 || roadsRet[2*i] >= NC + NJ ||
                    roadsRet[2*i+1] < 0 || roadsRet[2*i+1] >= NC + NJ) {
                    addFatalError("You can only build roads between pairs of existing cities/junctions.");
                    return -1.0;
                }
                if (roadsRet[2*i] == roadsRet[2*i+1]) {
                    addFatalError("You can only build roads between distinct points.");
                    return -1.0;
                }
                if (roadsRet[2*i] >= NC && junctOk[roadsRet[2*i] - NC] == 0 || 
                    roadsRet[2*i+1] >= NC && junctOk[roadsRet[2*i+1] - NC] == 0) {
                    addFatalError("You can not build a road to a dysfunctional junction.");
                    return -1.0;
                }
                roadS[i] = Math.min(roadsRet[2*i], roadsRet[2*i+1]);
                roadE[i] = Math.max(roadsRet[2*i], roadsRet[2*i+1]);
                Integer road = new Integer(roadToInt(roadS[i], roadE[i]));
                if (rUsed.contains(road)) {
                    addFatalError("You can only build one road between each pair of points.");
                    return -1.0;
                }
                rUsed.add(road);
            }
        } else {
            // build no roads and no junctions
            NJ = NR = 0;
            junct = new Pnt[NJ];
            roadS = new int[NR];
            roadE = new int[NR];
        }

        // create a structure to handle cities and junctions uniformly
        vertices = new Pnt[NC + NJ];
        for (int i = 0; i < NC; ++i)
            vertices[i] = cities[i];
        for (int i = 0; i < NJ; ++i)
            vertices[NC + i] = junct[i];

        if (vis) {
            // draw the results of program execution (even with invalid score)
            draw();
        }

        // check that all cities and all functional junctions are connected by the roads
        java.util.List<ArrayList<Integer>> adj = new ArrayList<ArrayList<Integer>>(NC + NJ);
        for (int i = 0; i < NC + NJ; ++i)
            adj.add(new ArrayList<Integer>());
        for (int i = 0; i < NR; ++i) {
            adj.get(roadS[i]).add(roadE[i]);
            adj.get(roadE[i]).add(roadS[i]);
        }
        boolean[] used = new boolean[NC + NJ];
        ArrayList<Integer> next = new ArrayList<Integer>();
        next.add(0);
        while (!next.isEmpty()) {
            int vertexInd = next.get(0);
            next.remove(0);
            for (int i = 0; i < adj.get(vertexInd).size(); ++i) {
                int adjInd = adj.get(vertexInd).get(i);
                if (!used[adjInd]) {
                    used[adjInd] = true;
                    next.add(adjInd);
                }
            }
        }
        for (int i = 0; i < NC + NJ; ++i)
            if (!used[i] && (i < NC || i >= NC && junctOk[i - NC] == 1)) {
                addFatalError((i < NC ? ("City " + i) : ("Junction " + (i - NC)) ) + " not connected to the rest of the network.");
                return -1.0;
            }

        // score is a function of total length of the roads
        // and the number of junctions built
        double totalRoadLength = 0;
        for (int i = 0; i < NR; ++i)
            totalRoadLength += Math.sqrt(vertices[roadS[i]].dist2(vertices[roadE[i]]));
        if (debug) {
            addFatalError("Number of junctions built: " + NJ);
            addFatalError("Total length of roads built: " + totalRoadLength);
        }
        return totalRoadLength + jCost * NJ;
      }
      catch (Exception e) { 
        addFatalError("An exception occurred while trying to get your program's results.");
        e.printStackTrace(); 
        return -1.0;
      }
    }
// ------------- visualization part ------------
    static String exec, fileName;
    static boolean vis, debug;
    static Process proc;
    InputStream is;
    OutputStream os;
    BufferedReader br;
    // -----------------------------------------
    int[] buildJunctions(int S, int[] cities, double junctionCost, double failProb) throws IOException {
        StringBuffer sb = new StringBuffer();
        sb.append(S).append("\n");
        sb.append(cities.length).append("\n");
        for (int i = 0; i < cities.length; ++i) {
            sb.append(cities[i]).append("\n");
        }
        sb.append(junctionCost).append("\n");
        sb.append(failProb).append("\n");
        os.write(sb.toString().getBytes());
        os.flush();

        // and get the return value
        int N = Integer.parseInt(br.readLine());
        int[] ret = new int[N];
        for (int i = 0; i < N; i++)
            ret[i] = Integer.parseInt(br.readLine());
        return ret;
    }
    // -----------------------------------------
    int[] buildRoads(int[] junctionBuilt) throws IOException {
        StringBuffer sb = new StringBuffer();
        sb.append(junctionBuilt.length).append("\n");
        for (int i = 0; i < junctionBuilt.length; ++i) {
            sb.append(junctionBuilt[i]).append("\n");
        }
        os.write(sb.toString().getBytes());
        os.flush();

        // and get the return value
        int N = Integer.parseInt(br.readLine());
        int[] ret = new int[N];
        for (int i = 0; i < N; i++)
            ret[i] = Integer.parseInt(br.readLine());
        return ret;
    }
    // -----------------------------------------
    public void draw() {
        int SZX = S + 1, SZY = SZX;
        BufferedImage bi = new BufferedImage(SZX, SZY,BufferedImage.TYPE_INT_RGB);
        Graphics2D g2 = (Graphics2D)bi.getGraphics();

        // background
        g2.setColor(Color.WHITE);
        g2.fillRect(0, 0, SZX, SZY);

        // roads
        g2.setColor(Color.BLACK);
        for (int i = 0; i < NR; ++i) {
            g2.drawLine(vertices[roadS[i]].x, vertices[roadS[i]].y, vertices[roadE[i]].x, vertices[roadE[i]].y);
        }

        // cities and junctions as colored points
        g2.setColor(Color.BLACK);
        for (int i = 0; i < NC; ++i) {
            g2.fillOval(cities[i].x - 3, cities[i].y - 3, 7, 7);
        }
        for (int i = 0; i < NJ; ++i) {
            g2.setColor((junctOk[i] == 0 ? Color.RED : Color.GREEN));
            g2.fillOval(junct[i].x - 2, junct[i].y - 2, 5, 5);
        }

        try {
            ImageIO.write(bi,"png", new File(fileName + ".png"));
        } catch (Exception e) { e.printStackTrace(); }
    }
    // -----------------------------------------
    public RoadsAndJunctionsVis(String seed) {
      try {
        if (exec != null) {
            try {
                Runtime rt = Runtime.getRuntime();
                proc = rt.exec(exec);
                os = proc.getOutputStream();
                is = proc.getInputStream();
                br = new BufferedReader(new InputStreamReader(is));
                new ErrorReader(proc.getErrorStream()).start();
            } catch (Exception e) { e.printStackTrace(); }
        }
        System.out.println("Score = " + runTest(seed));
        if (proc != null)
            try { proc.destroy(); } 
            catch (Exception e) { e.printStackTrace(); }
      }
      catch (Exception e) { e.printStackTrace(); }
    }
    // -----------------------------------------
    public static void main(String[] args) {
        String seed = "1";
        vis = true;
        for (int i = 0; i<args.length; i++)
        {   if (args[i].equals("-seed"))
                seed = args[++i];
            if (args[i].equals("-exec"))
                exec = args[++i];
            if (args[i].equals("-novis"))
                vis = false;
            if (args[i].equals("-debug"))
                debug = true;
        }
        if (exec == null)
            vis = true;
        if (vis)
            fileName = seed;
        RoadsAndJunctionsVis f = new RoadsAndJunctionsVis(seed);
    }
    // -----------------------------------------
    void addFatalError(String message) {
        System.out.println(message);
    }
}

class ErrorReader extends Thread{
    InputStream error;
    public ErrorReader(InputStream is) {
        error = is;
    }
    public void run() {
        try {
            byte[] ch = new byte[50000];
            int read;
            while ((read = error.read(ch)) > 0)
            {   String s = new String(ch,0,read);
                System.out.print(s);
                System.out.flush();
            }
        } catch(Exception e) { }
    }
}
