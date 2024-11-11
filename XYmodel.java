import java.awt.*;
import java.awt.event.*;
import java.text.DecimalFormat;

class XYmodel extends Canvas implements Runnable {

    int size = 200;                             // number of lattice sites in a row (change if desired)
    int squareWidth = 5;                        // pixels across one lattice site (increased for better visibility)
    int canvasSize = size * squareWidth;        // total pixels across canvas
    double[][] s = new double[size][size];      // the 2D array of spins (angles between 0 and 2π)
    boolean running = false;                    // true when simulation is running
    boolean showArrows = false;                 // flag to toggle between color and arrow visualization
    Button startButton = new Button("  Start  ");
    Scrollbar tScroller;                        // scrollbar to adjust temperature
    Label tLabel = new Label("Temperature = 1.00  ");    // text label next to scrollbar
    DecimalFormat twoPlaces = new DecimalFormat("0.00");    // to format temperature readout
    Image offScreenImage;                       // for double-buffering
    Graphics offScreenGraphics;
    Checkbox arrowCheckbox;                     // checkbox to toggle visualization mode

    // Constructor method handles all the initializations:
    XYmodel() {
        Frame frame = new Frame("XY Model");       // initialize the GUI...
        frame.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                System.exit(0);                            // close button exits program
            }
        });
        Panel canvasPanel = new Panel();
        frame.add(canvasPanel);
        canvasPanel.add(this);
        setSize(canvasSize,canvasSize);
        Panel controlPanel = new Panel();
        frame.add(controlPanel,BorderLayout.SOUTH);
        controlPanel.add(tLabel);
        tScroller = new Scrollbar(Scrollbar.HORIZONTAL,100,1,1,1001) {
            public Dimension getPreferredSize() {
                return new Dimension(100,15);            // make it bigger than default
            }
        };
        tScroller.setBlockIncrement(1);        // enables fine adjustments
        tScroller.addAdjustmentListener(new AdjustmentListener() {
            public void adjustmentValueChanged(AdjustmentEvent e) {
                tLabel.setText("Temperature = " + twoPlaces.format(tScroller.getValue()/100.0));
            }
        });
        controlPanel.add(tScroller);
        controlPanel.add(new Label("     "));            // leave some space

        // Add the checkbox to toggle visualization mode
        arrowCheckbox = new Checkbox("Show Arrows");
        arrowCheckbox.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                showArrows = arrowCheckbox.getState();
                // Repaint the entire lattice to update the visualization
                for (int i=0; i < size; i++) {
                    for (int j=0; j < size; j++) {
                        colorSquare(i,j);
                    }
                }
                repaint();
            }
        });
        controlPanel.add(arrowCheckbox);
        controlPanel.add(new Label("     "));            // leave some space

        startButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                running = !running;
                if (running) startButton.setLabel("Pause"); else startButton.setLabel("Resume");
            }
        });
        controlPanel.add(startButton);
        frame.pack();
        offScreenImage = createImage(canvasSize,canvasSize);
        offScreenGraphics = offScreenImage.getGraphics();

        for (int i=0; i < size; i++) {                    // initialize the lattice...
            for (int j=0; j < size; j++) {
                s[i][j] = Math.random() * 2 * Math.PI;    // random angle between 0 and 2π
                colorSquare(i,j);
            }
        }

        frame.setVisible(true);          // we're finally ready to show it!

        Thread t = new Thread(this);          // create a thread to run the simulation
        t.start();                            // and let 'er rip...
    }

    // Run method gets called by new thread to carry out the simulation:
    public void run() {
        double deltaTheta = Math.PI;  // maximum change in angle
        while (true) {
            if (running) {
                double temp = tScroller.getValue() / 100.0;
                for (int step=0; step<10000; step++) {       // adjust number of steps as desired
                    int i = (int) (Math.random() * size);    // choose a random row and column
                    int j = (int) (Math.random() * size);
                    double thetaOld = s[i][j];
                    double deltaThetaRandom = (Math.random() - 0.5) * 2 * deltaTheta; // random change in [-deltaTheta, deltaTheta]
                    double thetaNew = thetaOld + deltaThetaRandom;
                    thetaNew = thetaNew % (2 * Math.PI);
                    if (thetaNew < 0) thetaNew += 2 * Math.PI;

                    double eDiff = deltaU(i, j, thetaNew);

                    if ((eDiff <= 0) || (Math.random() < Math.exp(-eDiff / temp))) {    // Metropolis!
                        s[i][j] = thetaNew;
                        colorSquare(i,j);
                    }
                }
                repaint();        // causes update method to be called soon
            }
            try { Thread.sleep(1); } catch (InterruptedException e) {}  // sleep time in milliseconds
        }
    }

    // Given a lattice site and proposed new angle, compute energy change; note pbc:
    double deltaU(int i, int j, double thetaNew) {
        double thetaOld = s[i][j];
        double leftTheta, rightTheta, topTheta, bottomTheta;

        if (i == 0) leftTheta = s[size-1][j]; else leftTheta = s[i-1][j];
        if (i == size-1) rightTheta = s[0][j]; else rightTheta = s[i+1][j];
        if (j == 0) topTheta = s[i][size-1]; else topTheta = s[i][j-1];
        if (j == size-1) bottomTheta = s[i][0]; else bottomTheta = s[i][j+1];

        double deltaU = - (Math.cos(thetaNew - leftTheta) + Math.cos(thetaNew - rightTheta) +
                           Math.cos(thetaNew - topTheta) + Math.cos(thetaNew - bottomTheta))
                        + (Math.cos(thetaOld - leftTheta) + Math.cos(thetaOld - rightTheta) +
                           Math.cos(thetaOld - topTheta) + Math.cos(thetaOld - bottomTheta));
        return deltaU;
    }

    // Color a given square or draw an arrow according to the site's orientation:
    void colorSquare(int i, int j) {
        if (showArrows) {
            // Draw a white background
            offScreenGraphics.setColor(Color.white);
            offScreenGraphics.fillRect(i * squareWidth, j * squareWidth, squareWidth, squareWidth);

            // Draw an arrow representing the spin direction
            double theta = s[i][j] % (2 * Math.PI);

            int xCenter = i * squareWidth + squareWidth / 2;
            int yCenter = j * squareWidth + squareWidth / 2;
            int arrowLength = (int) (squareWidth * 0.8); // arrow length is 80% of square size
            int halfLength = arrowLength / 2;

            int xStart = (int) (xCenter - halfLength * Math.cos(theta));
            int yStart = (int) (yCenter + halfLength * Math.sin(theta));
            int xEnd = (int) (xCenter + halfLength * Math.cos(theta));
            int yEnd = (int) (yCenter - halfLength * Math.sin(theta));

            // Draw the arrow
            if (theta < 0)
                theta += 2 * Math.PI;
            float hue = (float) (theta / (2 * Math.PI));
            Color c = Color.getHSBColor(hue, 1.0f, 1.0f);
            offScreenGraphics.setColor(c);
            drawArrow(offScreenGraphics, xStart, yStart, xEnd, yEnd, squareWidth / 5, squareWidth / 5);

        } else {
            // Use the colored square visualization
            double theta = s[i][j] % (2 * Math.PI);
            if (theta < 0) theta += 2 * Math.PI;
            float hue = (float)(theta / (2 * Math.PI));
            Color c = Color.getHSBColor(hue, 1.0f, 1.0f);
            offScreenGraphics.setColor(c);
            offScreenGraphics.fillRect(i*squareWidth, j*squareWidth, squareWidth, squareWidth);
        }
    }

    // Method to draw an arrow from (x1, y1) to (x2, y2)
    void drawArrow(Graphics g, int x1, int y1, int x2, int y2, int d, int h) {
        int dx = x2 - x1, dy = y2 - y1;
        double D = Math.sqrt(dx*dx + dy*dy);
        double xm = D - d, xn = xm, ym = h, yn = -h, x;
        double sin = dy / D, cos = dx / D;
    
        x = xm*cos - ym*sin + x1;
        ym = xm*sin + ym*cos + y1;
        xm = x;
    
        x = xn*cos - yn*sin + x1;
        yn = xn*sin + yn*cos + y1;
        xn = x;
    
        int[] xpoints = {x2, (int) xm, (int) xn};
        int[] ypoints = {y2, (int) ym, (int) yn};
    
        g.drawLine(x1, y1, x2, y2);
        g.fillPolygon(xpoints, ypoints, 3);
    }

    // Override default update method to skip drawing the background:
    public void update(Graphics g) {
        paint(g);
    }

    // Paint method just blasts the off-screen image to the screen:
    public void paint(Graphics g) {
        g.drawImage(offScreenImage,0,0,canvasSize,canvasSize,this);
    }

    // Main method just calls constructor to get started:
    public static void main(String[] args) {
        new XYmodel();    
    }  
}
