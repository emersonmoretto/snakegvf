import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FileDialog;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.color.ColorSpace;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.DecimalFormat;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.WindowConstants;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class SnakeGUI {

	// --- MODEL -----------------------------------------------------------

	// the true snake object
	private Snake snakeinstance = null;

	// ---- IMAGE DATA -----------------------------------------------------

	private BufferedImage image = null;
	private BufferedImage imageanimation = null;
	private int[][] chanel_gradient = null;
	private double[][] chanel_flow = null;

	private int[][] intensite_flux = null;
	private int[][] r = null;
	// --- SWING COMPONENTS ------------------------------------------------

	private JLabel label0 = new JLabel("Gradient Vector Flow");
	private JLabel label1 = new JLabel("Snake");
	private JCheckBox cbShowAnim = new JCheckBox("Show animation");
	private JTextField txtMaxiter = new JTextField("500", 3);

	private JSlider slideThreshold = new JSlider(1, 100, 25);
	private JTextField txtAlpha = new JTextField("1.0", 3);
	private JTextField txtBeta = new JTextField("0.2", 3);
	private JTextField txtGamma = new JTextField("1.0", 3);
	private JTextField txtDelta = new JTextField("1.0", 3);

	private JCheckBox cbAutoadapt = new JCheckBox("Auto-Adapt");
	private JTextField txtStep = new JTextField("10", 3);
	private JTextField txtMinlen = new JTextField("2", 3);
	private JTextField txtMaxlen = new JTextField("3", 3);

	// --- Create the GUI ---------------------------------------------------

	public SnakeGUI() {
		label0.setVerticalTextPosition(JLabel.BOTTOM);
		label0.setHorizontalTextPosition(JLabel.CENTER);
		label1.setVerticalTextPosition(JLabel.BOTTOM);
		label1.setHorizontalTextPosition(JLabel.CENTER);

		final JFrame frame = new JFrame("Snake GVF - emoretto at usp br");
		frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);

		final JPanel panneau = new JPanel();
		panneau.add(label0);
		panneau.add(label1);
		final JScrollPane scrollPane = new JScrollPane(panneau);

		JButton buttonLoad = new JButton("Load");
		JButton buttonRun = new JButton("Run");

		final JPanel buttonPanel = new JPanel();
		buttonPanel.add(buttonLoad);
		buttonPanel.add(buttonRun);
		buttonPanel.add(cbShowAnim);
		buttonPanel.add(new JLabel("Max. Iteration:"));
		buttonPanel.add(txtMaxiter);
		cbShowAnim.setSelected(true);

		final JPanel coefPanel = new JPanel();
		coefPanel.add(new JLabel("gradient threshold:"));
		coefPanel.add(slideThreshold);
		coefPanel.add(new JLabel("alpha:"));
		coefPanel.add(txtAlpha);
		coefPanel.add(new JLabel("beta:"));
		coefPanel.add(txtBeta);
		coefPanel.add(new JLabel("gamma:"));
		coefPanel.add(txtGamma);
		coefPanel.add(new JLabel("delta:"));
		coefPanel.add(txtDelta);

		final JPanel adpatPanel = new JPanel();
		adpatPanel.add(cbAutoadapt);
		adpatPanel.add(new JLabel("every X iterations:"));
		adpatPanel.add(txtStep);
		adpatPanel.add(new JLabel("min segment length:"));
		adpatPanel.add(txtMinlen);
		adpatPanel.add(new JLabel("max segment length:"));
		adpatPanel.add(txtMaxlen);
		cbAutoadapt.setSelected(true);

		final JPanel mainPanel = new JPanel();
		mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.Y_AXIS));
		mainPanel.add(buttonPanel);
		mainPanel.add(coefPanel);
		mainPanel.add(adpatPanel);

		// frame
		frame.getContentPane().add(scrollPane, BorderLayout.CENTER);
		frame.getContentPane().add(mainPanel, BorderLayout.PAGE_END);
		frame.setSize(700, 500);
		frame.setVisible(true);

		// ActionListener "LOAD"
		buttonLoad.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				FileDialog filedialog = new FileDialog(frame, "Choose an image file");
				filedialog.setVisible(true);
				String filename = filedialog.getFile();
				String directory = filedialog.getDirectory();
				if (filename == null)
					return;
				File file = new File(directory + File.separator + filename);
				mainPanel.setVisible(false);
				try {
					loadimage(file);
					computegflow();
				} catch (Exception ex) {
					error(ex.getMessage(), ex);
				}
				mainPanel.setVisible(true);
			}
		});

		// ActionListener "SLIDE"
		slideThreshold.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if (slideThreshold.getValueIsAdjusting())
					return;
				mainPanel.setVisible(false);
				try {
					computegflow();
				} catch (Exception ex) {
					error(ex.getMessage(), ex);
				}
				mainPanel.setVisible(true);
			}
		});

		// Snake Runnable
		final Runnable snakerunner = new Runnable() {
			public void run() {
				try {
					startsnake();
				} catch (Exception ex) {
					error(ex.getMessage(), ex);
				}
			}
		};

		// ActionListener "RUN"
		buttonRun.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				new Thread(snakerunner).start();
			}
		});

		/*
		try {
			mainPanel.setVisible(false);
			loadimage(new File("trefle.png"));
			computegflow();
			mainPanel.setVisible(true);
		} catch (Exception ex) {
			error(ex.getMessage(), ex);
		}*/
	}

	// error (exception) display
	private static void error(String text, Exception ex) {
		if (ex != null) {
			StringWriter sw = new StringWriter();
			PrintWriter pw = new PrintWriter(sw);
			ex.printStackTrace(pw);
			String s = sw.toString();
			s = s.substring(0, Math.min(512, s.length()));
			text = text + "\n\n" + s + " (...)";
		}
		JOptionPane.showMessageDialog(null, text, "Error", JOptionPane.ERROR_MESSAGE);
	}
	
	// ---------------------------------------------------------------------
	//                        DRAWING PRIMITIVES
	// ---------------------------------------------------------------------

	public void display() { /* callback from snakeinstance */
		Graphics2D gc = imageanimation.createGraphics();
		
		// draw background image
        gc.drawImage(image,0,0,null);

		// draw snake lines
		gc.setColor( Color.RED );
		gc.setStroke(new BasicStroke(2.0f));
		List<Point> snakepoints = snakeinstance.snake;
		for (int i = 0; i < snakepoints.size(); i++) {
			int j = (i + 1) % snakepoints.size();
			Point p1 = snakepoints.get(i);
			Point p2 = snakepoints.get(j);
			gc.drawLine(p1.x, p1.y, p2.x, p2.y);
		}

		// draw snake points
		gc.setColor( Color.ORANGE );
		for (int i = 0; i < snakepoints.size(); i++) {
			Point p = snakepoints.get(i);
			gc.fillRect(p.x-2, p.y-2, 5, 5);
		}

		// swing display
		label1.setIcon(new ImageIcon(imageanimation));
	}

	// ---------------------------------------------------------------------
	//                       IMAGE LOADING/COMPUTATION
	// ---------------------------------------------------------------------

	private void loadimage(File file) throws IOException {
		image = null;
		image = ImageIO.read( file );
		imageanimation = new BufferedImage(image.getWidth(),image.getHeight(),ColorSpace.TYPE_RGB);

		// swing display
		label1.setIcon(new ImageIcon(image));
	}

	 
	/**
	 * @param f : image normalized in [0,1] 
	 * @param w : width of image
	 * @param h : height of image
	 * @param ITER : number of iterations
	 * @param mu : iteration step
	 * @return u[x,y] and v[x,y] arrays
	 */
	public double[][][] gvf(double[][] f, int w, int h, int ITER, double mu) {
	 
		// create empty arrays
		double[][] u = new double[w][h];
		double[][] v = new double[w][h];
		double[][] fx = new double[w][h];
		double[][] fy = new double[w][h];
		double[][] Lu = new double[w][h];
		double[][] Lv = new double[w][h];
	 
		// precompute edge-map (gradient)
		for (int y=1;y<(h-1);y++) {
			for (int x=1;x<(w-1);x++) {
				fx[x][y] = (f[x+1][y]-f[x-1][y])/2;
				fy[x][y] = (f[x][y+1]-f[x][y-1])/2;
				u[x][y] = fx[x][y];
				v[x][y] = fy[x][y];
			}
		}
		
	 
		// iterative diffusion
		for(int loop=0;loop<ITER;loop++) {
	 
			// compute laplacian of U and V
			for (int y=0 ; y < h ; y++) {
				for (int x=0 ; x < w ; x++) {
					
					if(x > 0 && y > 0 && x < w-1 && y < h-1){
						//Lu[x][y] = 0.25*(u[x-1][y]+u[x+1][y]+u[x][y-1]+u[x][y+1]) -u[x][y]; 
						//Lv[x][y] = 0.25*(v[x-1][y]+v[x+1][y]+v[x][y-1]+v[x][y+1]) -v[x][y];
						
						Lu[x][y] = ((u[x-1][y]+u[x+1][y] + u[x][y-1] + u[x][y+1]) - 4 * u[x][y]) / 4; 
						Lv[x][y] = ((v[x-1][y]+v[x+1][y] + v[x][y-1] + v[x][y+1]) - 4 * v[x][y]) / 4;
						
						//interior(i,j) = ((in(i-1,j) + in(i+1,j) + in(i,j-1) + in(i,j+1)) - 4*in(i,j) )/(h**2)
							    
								
					}else{

						if(x==0 && y ==0){
							Lu[x][y] = -5 * u[x][y+1] + 4 * u[x][y+2] - u[x][y+3] + 2 * u[x][y] - 5 * u[x+1][y] + 4 * u[x+2][y] - u[x+3][y] + 2 * u[x][y] / 16;
							Lv[x][y] = -5 * v[x][y+1] + 4 * v[x][y+2] - v[x][y+3] + 2 * v[x][y] - 5 * v[x+1][y] + 4 * v[x+2][y] - v[x+3][y] + 2 * v[x][y] / 16;
							Lu[x][y] /= 4;
							Lv[x][y] /= 4;
						}
					
						else if(y == 0 && x < w-1){
							Lu[x][y] = -5 * u[x][y+1] + 4 * u[x][y+2] - u[x][y+3] + 2 * u[x][y] + u[x+1][y] + u[x-1][y] - 2 * u[x][y] / 16;
							Lv[x][y] = -5 * v[x][y+1] + 4 * v[x][y+2] - v[x][y+3] + 2 * v[x][y] + v[x+1][y] + v[x-1][y] - 2 * v[x][y] / 16;
							Lu[x][y] /= 4;
							Lv[x][y] /= 4;
						}
						
						else if(x == 0  && y < h-1){
							Lu[x][y] = -5 * u[x+1][y] + 4 * u[x+2][y] - u[x+3][y] + 2 * u[x][y] + u[x][y+1] + u[x][y-1] - 2 * u[x][y] / 16;
							Lv[x][y] = -5 * v[x+1][y] + 4 * v[x+2][y] - v[x+3][y] + 2 * v[x][y] + v[x][y+1] + v[x][y-1] - 2 * v[x][y] / 16;
							Lu[x][y] /= 4;
							Lv[x][y] /= 4;
						}
						else if(x == w-1 && y > 0 && y < h-1){
							Lu[x][y] = -5 * u[x-1][y] + 4 * u[x-2][y] - u[x-3][y] + 2 * u[x][y] + u[x][y+1] + u[x][y-1] - 2 * u[x][y] / 16;
							Lv[x][y] = -5 * v[x-1][y] + 4 * v[x-2][y] - v[x-3][y] + 2 * v[x][y] + v[x][y+1] + v[x][y-1] - 2 * v[x][y] / 16;
							Lu[x][y] /= 4;
							Lv[x][y] /= 4;
						}
						else if(y == h-1 && x > 0  && x < w-1){
							Lu[x][y] = -5 * u[x][y-1] + 4 * u[x][y-2] - u[x][y-3] + 2 * u[x][y] + u[x+1][y] + u[x-1][y] - 2 * u[x][y] / 16;
							Lv[x][y] = -5 * v[x][y-1] + 4 * v[x][y-2] - v[x][y-3] + 2 * v[x][y] + v[x+1][y] + v[x-1][y] - 2 * v[x][y] / 16;
							Lu[x][y] /= 4;
							Lv[x][y] /= 4;
						}
						
					}
					
					//Lu[x][y] /= 4;
					//Lv[x][y] /= 4;
				}
			}
			
			debug(Lv);
//			try {
//				System.in.read();
//			} catch (IOException e) {
//			}
			
			
			// Calculate BC for the corners
			/*
			double ul = (-5.0*Lu[1][2] + 4.0*in(1,3) - in(1,4) + 2.0*in(1,1) - 5.0*in(2,1) + 4.0*in(3,1) - in(4,1) + 2.0*in(1][1))/4;
			double br = (-5.0*in(nr][nc-1) + 4.0*in(nr,nc-2) - in(nr,nc-3) + 2.0*in(nr,nc) - 5.0*in(nr-1,nc) + 4.0*in(nr-2,nc) - in(nr-3][nc) + 2.0*in(nr,nc))/4;
			double bl = (-5.0*in(nr][2) + 4.0*in(nr,3) - in(nr,4) + 2.0*in(nr,1) - 5.0*in(nr-1,1) + 4.0*in(nr-2,1) - in(nr-3,1) + 2.0*in(nr][1))/4;
			double ur = (-5.0*in(1][nc-1) + 4.0*in(1,nc-2) - in(1,nc-3) + 2.0*in(1,nc) - 5.0*in(2,nc) + 4.0*in(3,nc) - in(4,nc) + 2.0*in(1][nc))/4;
			*/
			
			
			/*
			mag = sqrt(u.*u+v.*v);
			px = u./(mag+1e-10); 
			py = v./(mag+1e-10);
			*/
			// update U and V
			for (int y=0;y<h;y++) {
				for (int x=0;x<w;x++) {
					double gnorm2 = fx[x][y]*fx[x][y] + fy[x][y]*fy[x][y];
	 
					u[x][y] += mu*4*Lu[x][y] - gnorm2 * (u[x][y]-fx[x][y]);
					v[x][y] += mu*4*Lv[x][y] - gnorm2 * (v[x][y]-fy[x][y]);
					
					chanel_flow[x][y] = Math.sqrt(u[x][y]*u[x][y] + v[x][y]*v[x][y]);
					
					//chanel_flow[x][y] = (int) u[y][y] + (int) v[x][y] ;
					//chanel_flow[x][y] = u[x][y] + v[x][y] / ( Math.sqrt(u[x][y]*u[x][y] + v[x][y]*v[x][y]) + 1e-10 );
					//chanel_flow[x][y] += v[x][y] / Math.sqrt(u[x][y]*u[x][y] + v[x][y]*v[x][y]) + 1e-10;
				}
			}
			
			
		}
		
		
	 
		// return U and V arrays
		return new double[][][]{u,v};
	}
	
	private void normalize2(double[][] mtx){
			
		/*
		mag = sqrt(u.*u+v.*v);
		px = u./(mag+1e-10); 
		py = v./(mag+1e-10);
		*/
		
	//101	4
		
	}
	
	private String debug(double[][] mtx) {
	 
		DecimalFormat df = new DecimalFormat("#.####");
		StringBuilder sb = new StringBuilder();
		for(int i=0 ; i < mtx.length ; i++){
			for(int j=0 ; j < mtx[i].length ; j++){
				System.out.print(df.format(mtx[i][j])+ "\t");
				sb.append(df.format(mtx[i][j])+ "\t");
				//mtx[i][j] = mtx[i][j] / Math.sqrt() + a;
				//mtx[i][j] /= sum;
			}
			sb.append("\n");
			System.out.println("");
		}
		return sb.toString();	
			
	}

	// normalize energy matrix
	private void normalize(double[][] mtx) {
		
		double a = 1e-10;
		double sum=0;		
		for(int i=0 ; i < mtx.length ; i++)
			for(int j=0 ; j < mtx[i].length ; j++)
				sum += Math.abs(mtx[i][j]);
	 
		if (sum==0) return;
	 
		for(int i=0 ; i < mtx.length ; i++)
			for(int j=0 ; j < mtx[i].length ; j++){
				
				
				//mtx[i][j] = mtx[i][j] / Math.sqrt() + a;
				//mtx[i][j] /= sum;
			}
	}

	
	private void computegflow() {
		int W = image.getWidth();
		int H = image.getHeight();

		this.intensite_flux = new int[W][H];
		this.r = new int[W][H]; 
		
		int THRESHOLD = slideThreshold.getValue();

		// GrayLevelScale (Luminance)
		int[][] clum = new int[W][H];
		for (int y = 0; y < H; y++)
			for (int x = 0; x < W; x++) {
				int rgb=image.getRGB(x,y);
				int r = (rgb >>16 ) & 0xFF;
				int g = (rgb >> 8 ) & 0xFF;
				int b = rgb & 0xFF;
				clum[x][y] = (int)(0.299*r + 0.587*g + 0.114*b);  
			}
				
		// Gradient (sobel)
		this.chanel_gradient = new int[W][H]; 
		int maxgradient=0;
		for (int y = 0; y < H-2; y++)
			for (int x = 0; x < W-2; x++) {
				int p00 = clum[x+0][y+0]; int p10 = clum[x+1][y+0]; int p20 = clum[x+2][y+0];
				int p01 = clum[x+0][y+1]; /*-------------------- */ int p21 = clum[x+2][y+1];
				int p02 = clum[x+0][y+2]; int p12 = clum[x+1][y+2]; int p22 = clum[x+2][y+2];
				int sx = (p20+2*p21+p22)-(p00+2*p01+p02);
				int sy = (p02+2*p12+p22)-(p00+2*p10+p10);
				int snorm = (int)Math.sqrt(sx*sx+sy*sy);
				chanel_gradient[x+1][y+1]=snorm;
				maxgradient=Math.max(maxgradient, snorm);
			}

		// thresholding
		boolean[][] binarygradient = new boolean[W][H];
		for (int y = 0; y < H; y++)
			for (int x = 0; x < W; x++)
				if (chanel_gradient[x][y] > THRESHOLD*maxgradient/100) {
					binarygradient[x][y]=true;
				} else {
					chanel_gradient[x][y]=0;
				}
		// distance map to binarized gradient
		chanel_flow = new double[W][H];
		
		double[][] cdist = new ChamferDistance(ChamferDistance.chamfer5).compute(binarygradient, W,H);
		for (int y = 0; y < H; y++)
			for (int x = 0; x < W; x++)
				chanel_flow[x][y]=(int)(5*cdist[x][y]);
		
		
		/**
		 * to gray and normalize
		 */
		double[][] f = new double[W][H];
		
		for(int i=0 ; i < image.getWidth() ; i++)
			for(int j=0 ; j < image.getHeight() ; j++){
				
				 int rgbP = image.getRGB(i, j);
			     int r = (rgbP >> 16) & 0xFF;
			     int g = (rgbP >> 8) & 0xFF;
			     int b = (rgbP & 0xFF);

			     int grayLevel = (r + g + b) / 3;
			        
			     //norm
			     f[i][j] = grayLevel / 255;
			}
		
		
		
				
		// THE GVF!!!
		double[][][] gvfield = gvf(f, W, H, 50, 0.2);
		
		debug(chanel_flow);
		/*
		//normalize(f);
		for (int y=0;y<H;y++) {
			for (int x=0;x<W;x++) {
				//chanel_flow[x][y] = Math.sqrt(u[x][y]*u[x][y] + v[x][y]*v[x][y]);
				//chanel_flow[x][y] = (255 - (chanel_flow[x][y]*1000) );
				System.out.print( ((double) chanel_flow[x][y] + " ") );
			}		
			System.out.println("");
		}*/
		
		//normalize(chanel_flow);
		
		
		// show flow + gradient
		int[] rgb = new int[3];
		BufferedImage imgflow = new BufferedImage(W, H, ColorSpace.TYPE_RGB);
		
		for (int y = 0; y < H; y++) {
			for (int x = 0; x < W; x++) {
				int vflow = (int) chanel_flow[x][y];
				int vgrad = binarygradient[x][y]?255:0;

				if (vgrad > 0) {
					rgb[0] = 0;
					rgb[1] = vgrad;
					rgb[2] = 0;
				} else {
					rgb[0] = Math.max(0, 255 - vflow);
					rgb[1] = 0;
					rgb[2] = 0;
				}
				int irgb = (0xFF<<24)+(rgb[0]<<16)+(rgb[1]<<8)+rgb[2];
				imgflow.setRGB(x, y, irgb);
			}
		}
		
		

		// swing display
		label0.setIcon(new ImageIcon(imgflow));
	}

	// ---------------------------------------------------------------------
	//                         START SNAKE SEGMENTATION
	// ---------------------------------------------------------------------

	private void startsnake() {
		int W = image.getWidth();
		int H = image.getHeight();
		int MAXLEN = Integer.parseInt(txtMaxlen.getText()); /* max segment length */

		// initial points
		double radius = (W/2 + H/2) / 2;
		double perimeter = 6.28 * radius;
		int nmb = (int) (perimeter / MAXLEN);
		Point[] circle = new Point[nmb];
		for (int i = 0; i < circle.length; i++) {
			double x = (W / 2 + 0) + (W / 2 - 2) * Math.cos((6.28 * i) / circle.length);
			double y = (H / 2 + 0) + (H / 2 - 2) * Math.sin((6.28 * i) / circle.length);
			circle[i] = new Point((int) x, (int) y);
		}

		// create snake instance
		snakeinstance = new Snake(W, H, chanel_gradient, chanel_flow, circle);

		// snake base parameters
		snakeinstance.alpha = Double.parseDouble(txtAlpha.getText());
		snakeinstance.beta = Double.parseDouble(txtBeta.getText());
		snakeinstance.gamma = Double.parseDouble(txtGamma.getText());
		snakeinstance.delta = Double.parseDouble(txtDelta.getText());

		// snake extra parameters
		snakeinstance.SNAKEGUI = this;
		snakeinstance.SHOWANIMATION = cbShowAnim.isSelected();
		snakeinstance.AUTOADAPT = cbAutoadapt.isSelected();
		snakeinstance.AUTOADAPT_LOOP = Integer.parseInt(txtStep.getText());
		snakeinstance.AUTOADAPT_MINLEN = Integer.parseInt(txtMinlen.getText());
		snakeinstance.AUTOADAPT_MAXLEN = Integer.parseInt(txtMaxlen.getText());
		snakeinstance.MAXITERATION = Integer.parseInt(txtMaxiter.getText());

		// animate snake
		System.out.println("initial snake points:" + snakeinstance.snake.size());
		int nmbloop = snakeinstance.loop();
		System.out.println("final snake points:" + snakeinstance.snake.size());
		System.out.println("iterations: " + nmbloop);

		// display final result
		display();

		System.out.println("END");
	}

	// ---------------------------------------------------------------------

	public static void main(String[] args) {
		Thread.currentThread().setPriority(Thread.MIN_PRIORITY);
		new SnakeGUI();
	}
}
