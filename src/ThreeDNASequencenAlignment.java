import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ThreeDNASequencenAlignment {
	private int[][] blast_scores;
	private String alphabet = "ACGT-";
	private String s1,s2,s3;
	private List<Integer> alignment_scores;
	private StringBuilder sb1, sb2, sb3;
	private String rs1, rs2, rs3;
	private int s1_lenth, s2_length, s3_length;
	
	public ThreeDNASequencenAlignment(String fn1, String fn2, String fn3) throws IOException
	{
		this.s1 = read_file(fn1);
		this.s2 = read_file(fn2);
		this.s3 = read_file(fn3);
		this.s1_lenth = s1.length();
		this.s2_length = s2.length();
		this.s3_length = s3.length();
		System.out.println("Lenght of the first sequence: " + this.s1_lenth);
		System.out.println("Lenght of the second sequence: " + this.s2_length);
		System.out.println("Lenght of the third sequence: " + this.s3_length);
		
		this.rs1 = new StringBuilder(s1).reverse().toString();
		this.rs2 = new StringBuilder(s2).reverse().toString();
		this.rs3 = new StringBuilder(s3).reverse().toString();
		
		this.initialize_blast_score();
	}

	public void align() {
		if(this.alignment_scores == null)
		{
			int m = this.s1.length();
			int n = this.s2.length();
			int l = this.s3.length();
			
			this.sb1 = new StringBuilder();
			this.sb2 = new StringBuilder();
			this.sb3 = new StringBuilder();
			this.alignment_scores = new ArrayList<Integer>();
	
			long start_time = System.currentTimeMillis();
			if (m == 0 || n == 0 || l == 0)
			{
				dp_align(this.s1, this.s2, this.s3, this.sb1, this.sb2, this.sb3);
			}
			else
			{
				recur_align(0, 0, 0, m, n, l);
			}
			long end_time = System.currentTimeMillis();
			System.out.println("Total running time (seconds): " + (end_time - start_time)/1000);
		}
	}

	private void recur_align(int s1_start, int s2_start, int s3_start, int s1_end, int s2_end, int s3_end) 
	{
		if ((s3_start + 1 >= s3_end) || (s1_start + 1 >= s1_end) || (s2_start + 1 >= s2_end)) {
			StringBuilder base_sb1 = new StringBuilder();
			StringBuilder base_sb2 = new StringBuilder();
			StringBuilder base_sb3 = new StringBuilder();
			int s = dp_align(s1.substring(s1_start, s1_end), s2.substring(s2_start, s2_end), s3.substring(s3_start, s3_end), base_sb1, base_sb2, base_sb3);
			this.alignment_scores.add(s);
			sb1.append(base_sb1);
			sb2.append(base_sb2);
			sb3.append(base_sb3);
		} else {
			int mid = (s3_start + s3_end) / 2;
			int[][] forward_scores = compute_forward_3D_scores(s1, s2, s3, s1_start, s2_start, s3_start, s1_end, s2_end, mid);
			int[][] backward_scores = compute_backward_3D_scores(s1_start, s2_start, mid, s1_end, s2_end, s3_end);

			int i_max = 0, j_max = 0;
			int max_sum = forward_scores[i_max][j_max] + backward_scores[i_max][j_max];
			for (int i = 0; i < forward_scores.length; i++) {
				for (int j = 0; j < forward_scores[i].length; j++) {
					if (forward_scores[i][j] + backward_scores[i][j] > max_sum) {
						max_sum = forward_scores[i][j] + backward_scores[i][j];
						i_max = i;
						j_max = j;
					}
				}
			}
			this.recur_align(s1_start, s2_start, s3_start, s1_start + i_max, s2_start + j_max, mid);
			this.recur_align(s1_start + i_max, s2_start + j_max, mid, s1_end, s2_end, s3_end);
		}
	}

	private int dp_align(String s1, String s2, String s3, StringBuilder sb1, StringBuilder sb2, StringBuilder sb3) {
		int m = s1.length();
		int n = s2.length();
		int l = s3.length();
		int[][][] cur_aligned_scores = new int[m + 1][n + 1][l + 1];
		int[][][] prev_pos = new int[m + 1][n + 1][l + 1];
		
		for (int i = 0; i < prev_pos.length; i++) {
			for (int j = 0; j < prev_pos[i].length; j++) {
				for (int k = 0; k < prev_pos[i][j].length; k++) {
					prev_pos[i][j][k] = -1;
				}
			}
		}
		cur_aligned_scores[0][0][0] = 0;
		for (int i = 1; i < m + 1; i++) {
			cur_aligned_scores[i][0][0] = cur_aligned_scores[i - 1][0][0] + 2 * blast_scores[4][alphabet.indexOf(s1.charAt(i - 1))];
			prev_pos[i][0][0] = 0;
		}
		for (int j = 1; j < n + 1; j++) {
			cur_aligned_scores[0][j][0] = cur_aligned_scores[0][j - 1][0] + 2 * blast_scores[4][alphabet.indexOf(s2.charAt(j - 1))];
			prev_pos[0][j][0] = 1;
		}
		for (int k = 1; k < l + 1; k++) {
			cur_aligned_scores[0][0][k] = cur_aligned_scores[0][0][k - 1] + 2 * blast_scores[4][alphabet.indexOf(s3.charAt(k - 1))];
			prev_pos[0][0][k] = 2;
		}

		for (int i = 1; i < m + 1; i++) {
			for (int j = 1; j < n + 1; j++) {
				int[] score = new int[3];
				List<Integer> scoreList = new ArrayList<Integer>();
				
				score[0] = cur_aligned_scores[i - 1][j - 1][0] + blast_scores[alphabet.indexOf(s1.charAt(i - 1))][alphabet.indexOf(s2.charAt(j - 1))]
							+ blast_scores[4][alphabet.indexOf(s1.charAt(i - 1))] 
							+ blast_scores[4][alphabet.indexOf(s2.charAt(j - 1))];
				score[1] = cur_aligned_scores[i - 1][j][0] + 2 * blast_scores[4][alphabet.indexOf(s1.charAt(i - 1))];
				score[2] = cur_aligned_scores[i][j - 1][0] + 2 * blast_scores[4][alphabet.indexOf(s2.charAt(j - 1))];
				
				scoreList.add(score[0]);
				scoreList.add(score[1]);
				scoreList.add(score[2]);
				
				int max_score = Collections.max(scoreList);
				int max_index = scoreList.indexOf(max_score);
				if (max_index == 0)
				{
					prev_pos[i][j][0] = 3;
				}
				else if (max_index == 1)
				{
					prev_pos[i][j][0] = 0;
				}
				else if (max_index == 2)
				{
					prev_pos[i][j][0] = 1;
				}
				cur_aligned_scores[i][j][0] = max_score;
			}
		}

		for (int i = 1; i < m + 1; i++) {
			for (int k = 1; k < l + 1; k++) {
				int[] score = new int[3];
				List<Integer> scoreList = new ArrayList<Integer>();
				
				score[0] = cur_aligned_scores[i - 1][0][k - 1] + blast_scores[alphabet.indexOf(s1.charAt(i - 1))][alphabet.indexOf(s3.charAt(k - 1))]
							+ blast_scores[4][alphabet.indexOf(s1.charAt(i - 1))] 
							+ blast_scores[4][alphabet.indexOf(s3.charAt(k - 1))];
				score[1] = cur_aligned_scores[i - 1][0][k] + 2 * blast_scores[4][alphabet.indexOf(s1.charAt(i - 1))];
				score[2] = cur_aligned_scores[i][0][k - 1] + 2 * blast_scores[4][alphabet.indexOf(s3.charAt(k - 1))];
				
				scoreList.add(score[0]);
				scoreList.add(score[1]);
				scoreList.add(score[2]);
				
				int max_score = Collections.max(scoreList);
				int max_index = scoreList.indexOf(max_score);
				if (max_index == 0)
				{
					prev_pos[i][0][k] = 4;
				}
				else if (max_index == 1)
				{
					prev_pos[i][0][k] = 0;
				}
				else if (max_index == 2)
				{
					prev_pos[i][0][k] = 2;
				}
				cur_aligned_scores[i][0][k] = max_score;
			}
		}

		for (int j = 1; j < n + 1; j++) {
			for (int k = 1; k < l + 1; k++) {
				int[] score = new int[3];
				List<Integer> scoreList = new ArrayList<Integer>();
				
				score[0] = cur_aligned_scores[0][j - 1][k - 1] + blast_scores[alphabet.indexOf(s2.charAt(j - 1))][alphabet.indexOf(s3.charAt(k - 1))]
							+ blast_scores[4][alphabet.indexOf(s2.charAt(j - 1))] 
							+ blast_scores[4][alphabet.indexOf(s3.charAt(k - 1))];
				score[1] = cur_aligned_scores[0][j][k - 1] + 2 * blast_scores[4][alphabet.indexOf(s3.charAt(k - 1))];
				score[2] = cur_aligned_scores[0][j - 1][k] + 2 * blast_scores[4][alphabet.indexOf(s2.charAt(j - 1))];
				
				scoreList.add(score[0]);
				scoreList.add(score[1]);
				scoreList.add(score[2]);
				
				int max_score = Collections.max(scoreList);
				int max_index = scoreList.indexOf(max_score);
				if (max_index == 0)
				{
					prev_pos[0][j][k] = 5;
				}
				else if (max_index == 1)
				{
					prev_pos[0][j][k] = 2;
				}
				else if (max_index == 2)
				{
					prev_pos[0][j][k] = 1;
				}
				cur_aligned_scores[0][j][k] = max_score;
			}
		}

		for (int k = 1; k < l + 1; k++) {
			for (int i = 1; i < m + 1; i++) {
				for (int j = 1; j < n + 1; j++) {
					int[] score = new int[7];
					List<Integer> scoreList = new ArrayList<Integer>();
					
					score[0] = cur_aligned_scores[i - 1][j][k] + 2 * blast_scores[4][alphabet.indexOf(s1.charAt(i - 1))];
					score[1] = cur_aligned_scores[i][j - 1][k] + 2 * blast_scores[4][alphabet.indexOf(s2.charAt(j - 1))];
					score[2] = cur_aligned_scores[i][j][k - 1] + 2 * blast_scores[4][alphabet.indexOf(s3.charAt(k - 1))];
					score[3] = cur_aligned_scores[i - 1][j - 1][k] + blast_scores[4][alphabet.indexOf(s1.charAt(i - 1))]
								+ blast_scores[4][alphabet.indexOf(s2.charAt(j - 1))]
								+ blast_scores[alphabet.indexOf(s1.charAt(i - 1))][alphabet.indexOf(s2.charAt(j - 1))];
					score[4] = cur_aligned_scores[i - 1][j][k - 1] + blast_scores[4][alphabet.indexOf(s1.charAt(i - 1))]
								+ blast_scores[4][alphabet.indexOf(s3.charAt(k - 1))]
								+ blast_scores[alphabet.indexOf(s1.charAt(i - 1))][alphabet.indexOf(s3.charAt(k - 1))];
					score[5] = cur_aligned_scores[i][j - 1][k - 1] + blast_scores[4][alphabet.indexOf(s2.charAt(j - 1))]
								+ blast_scores[4][alphabet.indexOf(s3.charAt(k - 1))]
								+ blast_scores[alphabet.indexOf(s2.charAt(j - 1))][alphabet.indexOf(s3.charAt(k - 1))];
					score[6] = cur_aligned_scores[i - 1][j - 1][k - 1] + blast_scores[alphabet.indexOf(s1.charAt(i - 1))][alphabet.indexOf(s2.charAt(j - 1))]
								+ blast_scores[alphabet.indexOf(s1.charAt(i - 1))][alphabet.indexOf(s3.charAt(k - 1))]
								+ blast_scores[alphabet.indexOf(s2.charAt(j - 1))][alphabet.indexOf(s3.charAt(k - 1))];
					
					scoreList.add(score[0]);
					scoreList.add(score[1]);
					scoreList.add(score[2]);
					scoreList.add(score[3]);
					scoreList.add(score[4]);
					scoreList.add(score[5]);
					scoreList.add(score[6]);
					int max_score = Collections.max(scoreList);
					int max_index = scoreList.indexOf(max_score);
					prev_pos[i][j][k] = max_index;
					cur_aligned_scores[i][j][k] = max_score;
					scoreList.clear();
				}
			}
		}

		trace_back(s1, s2, s3, prev_pos, m, n, l, sb1, sb2, sb3);
		sb1.reverse();
		sb2.reverse();
		sb3.reverse();
		return cur_aligned_scores[m][n][l];
	}

	private void trace_back(String s1, String s2, String s3, int[][][] prev_pos, int i, int j, int k,
			StringBuilder sb1, StringBuilder sb2, StringBuilder sb3) {
		if (i == 0 && j == 0 && k == 0) {
			return;
		}

		if (prev_pos[i][j][k] == 0) {
			sb1.append(s1.charAt(i - 1));
			sb2.append("-");
			sb3.append("-");
			trace_back(s1, s2, s3, prev_pos, i - 1, j, k, sb1, sb2, sb3);
		} else if (prev_pos[i][j][k] == 1) {
			sb1.append("-");
			sb2.append(s2.charAt(j - 1));
			sb3.append("-");
			trace_back(s1, s2, s3, prev_pos, i, j - 1, k, sb1, sb2, sb3);
		} else if (prev_pos[i][j][k] == 2) {
			sb1.append("-");
			sb2.append("-");
			sb3.append(s3.charAt(k - 1));
			trace_back(s1, s2, s3, prev_pos, i, j, k - 1, sb1, sb2, sb3);
		} else if (prev_pos[i][j][k] == 3) {
			sb1.append(s1.charAt(i - 1));
			sb2.append(s2.charAt(j - 1));
			sb3.append("-");
			trace_back(s1, s2, s3, prev_pos, i - 1, j - 1, k, sb1, sb2, sb3);
		} else if (prev_pos[i][j][k] == 4) {
			sb1.append(s1.charAt(i - 1));
			sb2.append("-");
			sb3.append(s3.charAt(k - 1));
			trace_back(s1, s2, s3, prev_pos, i - 1, j, k - 1, sb1, sb2, sb3);
		} else if (prev_pos[i][j][k] == 5) {
			sb1.append("-");
			sb2.append(s2.charAt(j - 1));
			sb3.append(s3.charAt(k - 1));
			trace_back(s1, s2, s3, prev_pos, i, j - 1, k - 1, sb1, sb2, sb3);
		} else if (prev_pos[i][j][k] == 6) {
			sb1.append(s1.charAt(i - 1));
			sb2.append(s2.charAt(j - 1));
			sb3.append(s3.charAt(k - 1));
			trace_back(s1, s2, s3, prev_pos, i - 1, j - 1, k - 1, sb1, sb2, sb3);
		} else
			System.out.println("Something wrong prev[i][j][k] = " + prev_pos[i][j][k]);
	}

	private int[][] compute_2D_scores(String s1, String s2, int s1_start, int s2_start, int s1_end, int s2_end) {
		int m = s1_end - s1_start;
		int n = s2_end - s2_start;
		int[][] scores = new int[m + 1][n + 1];

		scores[0][0] = 0;
		for (int i = 1; i < m + 1; i++) {
			scores[i][0] = scores[i-1][0] + 2 * blast_scores[4][alphabet.indexOf(s1.charAt(i - 1 + s1_start))];
		}

		for (int j = 1; j < n + 1; j++) {
			scores[0][j] = scores[0][j-1] + 2 * blast_scores[4][alphabet.indexOf(s2.charAt(j - 1 + s2_start))];
		}

		for (int i = 1; i < m + 1; i++) {
			for (int j = 1; j < n + 1; j++) {
				int[] score = new int[3];
				List<Integer> scoreList = new ArrayList<Integer>();
				
				score[0] = scores[i-1][j-1] + blast_scores[alphabet.indexOf(s1.charAt(i - 1 + s1_start))][alphabet.indexOf(s2.charAt(j - 1 + s2_start))]
							+ blast_scores[4][alphabet.indexOf(s1.charAt(i - 1 + s1_start))] 
							+ blast_scores[4][alphabet.indexOf(s2.charAt(j - 1 + s2_start))];
				score[1] = scores[i-1][j] + 2 * blast_scores[4][alphabet.indexOf(s1.charAt(i - 1 + s1_start))];
				score[2] = scores[i][j-1] + 2 * blast_scores[4][alphabet.indexOf(s2.charAt(j - 1 + s2_start))];
				
				scoreList.add(score[0]);
				scoreList.add(score[1]);
				scoreList.add(score[2]);
				int max_value = Collections.max(scoreList);
				scores[i][j] = max_value;
			}
		}
		return scores;
	}

	private int[][] compute_forward_3D_scores(String s1, String s2, String s3, int s1_start, int s2_start, int s3_start, int s1_end,
			int s2_end, int s3_end) {
		int[][] s12_scores = compute_2D_scores(s1, s2, s1_start, s2_start, s1_end, s2_end);
		int[][] s31_scores = compute_2D_scores(s3, s1, s3_start, s1_start, s3_end, s1_end);
		int[][] s32_scores = compute_2D_scores(s3, s2, s3_start, s2_start, s3_end, s2_end);

		int a, b, c, d, e, f, g;
		int[][] forward_scores = new int[s1_end - s1_start + 1][s2_end - s2_start + 1];
		for (int i = 0; i < s12_scores.length; i++) {
			for (int j = 0; j < s12_scores[i].length; j++) {
				forward_scores[i][j] = s12_scores[i][j];
			}
		}

		for (int k = 1; k < s3_end - s3_start + 1; k++) {
			/* Save scores of previous plane */
			for (int i = 0; i < s12_scores.length; i++) {
				for (int j = 0; j < s12_scores[i].length; j++) {
					s12_scores[i][j] = forward_scores[i][j];
				}
			}
			
			/* Initialize scores of current plane to all 0 */
			for (int i = 0; i < s12_scores.length; i++) {
				for (int j = 0; j < s12_scores[i].length; j++) {
					forward_scores[i][j] = 0;
				}
			}

			/* Assign initialization scores */
			for (int i = 0; i < s1_end - s1_start + 1; i++) {
				forward_scores[i][0] = s31_scores[k][i];
			}
			for (int j = 0; j < s2_end - s2_start + 1; j++) {
				forward_scores[0][j] = s32_scores[k][j];
			}

			for (int i = 1; i < s1_end - s1_start + 1; i++) {
				for (int j = 1; j < s2_end - s2_start + 1; j++) {
					int[] score = new int[7];
					List<Integer> scoreList = new ArrayList<Integer>();
					
					/* 4 points in previous plane */
					a = s12_scores[i - 1][j - 1];
					b = s12_scores[i - 1][j];
					c = s12_scores[i][j - 1];
					d = s12_scores[i][j];
					/* 3 points in current plane */
					e = forward_scores[i - 1][j - 1];
					f = forward_scores[i][j - 1];
					g = forward_scores[i - 1][j];
					
					score[0] = g + 2 * blast_scores[4][alphabet.indexOf(s1.charAt(i - 1 + s1_start))];
					score[1] = f + 2 * blast_scores[4][alphabet.indexOf(s2.charAt(j - 1 + s2_start))];
					score[2] = d + 2 * blast_scores[4][alphabet.indexOf(s3.charAt(k - 1 + s3_start))];
					score[3] = e + blast_scores[4][alphabet.indexOf(s1.charAt(i - 1 + s1_start))] 
								+ blast_scores[4][alphabet.indexOf(s2.charAt(j - 1 + s2_start))]
								+ blast_scores[alphabet.indexOf(s1.charAt(i - 1 + s1_start))][alphabet.indexOf(s2.charAt(j - 1 + s2_start))];
					score[4] = b + blast_scores[4][alphabet.indexOf(s1.charAt(i - 1 + s1_start))] 
								+ blast_scores[4][alphabet.indexOf(s3.charAt(k - 1 + s3_start))]
								+ blast_scores[alphabet.indexOf(s1.charAt(i - 1 + s1_start))][alphabet.indexOf(s3.charAt(k - 1 + s3_start))];
					score[5] = c + blast_scores[4][alphabet.indexOf(s2.charAt(j - 1 + s2_start))] 
								+ blast_scores[4][alphabet.indexOf(s3.charAt(k - 1 + s3_start))]
								+ blast_scores[alphabet.indexOf(s2.charAt(j - 1 + s2_start))][alphabet.indexOf(s3.charAt(k - 1 + s3_start))];
					score[6] = a + blast_scores[alphabet.indexOf(s1.charAt(i - 1 + s1_start))][alphabet.indexOf(s2.charAt(j - 1 + s2_start))]
								+ blast_scores[alphabet.indexOf(s1.charAt(i - 1 + s1_start))][alphabet.indexOf(s3.charAt(k - 1 + s3_start))]
								+ blast_scores[alphabet.indexOf(s2.charAt(j - 1 + s2_start))][alphabet.indexOf(s3.charAt(k - 1 + s3_start))];

					scoreList.add(score[0]);
					scoreList.add(score[1]);
					scoreList.add(score[2]);
					scoreList.add(score[3]);
					scoreList.add(score[4]);
					scoreList.add(score[5]);
					scoreList.add(score[6]);
					int max_value = Collections.max(scoreList);

					forward_scores[i][j] = max_value;
				}
			}
		}
		return forward_scores;
	}

	private int[][] compute_backward_3D_scores(int s1_start, int s2_start, int s3_start, int s1_end, int s2_end, int s3_end) {
		int rs1_start = this.s1_lenth - s1_end;
		int rs1_end = this.s1_lenth - s1_start;
		int rs2_start = this.s2_length - s2_end;
		int rs2_end = this.s2_length - s2_start;
		int rs3_start = this.s3_length - s3_end;
		int rs3_end = this.s3_length - s3_start;
		
		int[][] backward_scores = compute_forward_3D_scores(this.rs1, this.rs2, this.rs3, rs1_start, rs2_start, rs3_start, rs1_end, rs2_end, rs3_end);
		/* Rotate matrix 180 degree */
		rotate_matrix(backward_scores);

		return backward_scores;
	}

	private void rotate_matrix(int[][] matrix) {
		int start,end,temp;
		
		for (int i = 0; i < matrix.length; i++) {
			start = 0;
			end = matrix[i].length-1;
			while(start < end)
			{
				temp = matrix[i][start];
				matrix[i][start] = matrix[i][end];
				matrix[i][end] = temp;
				start++;
				end--;
			}
		}

		for (int j = 0; j < matrix[0].length; j++) {
			start = 0;
			end = matrix.length-1;
			while(start < end)
			{
				temp = matrix[start][j];
				matrix[start][j] = matrix[end][j];
				matrix[end][j] = temp;
				start++;
				end--;
			}
		}
	}

	private void initialize_blast_score() {
		blast_scores = new int[5][5];
		blast_scores[0][0] = 5;
		blast_scores[0][1] = -4;
		blast_scores[0][2] = -4;
		blast_scores[0][3] = -4;
		blast_scores[0][4] = -8;
		blast_scores[1][0] = -4;
		blast_scores[1][1] = 5;
		blast_scores[1][2] = -4;
		blast_scores[1][3] = -4;
		blast_scores[1][4] = -8;
		blast_scores[2][0] = -4;
		blast_scores[2][1] = -4;
		blast_scores[2][2] = 5;
		blast_scores[2][3] = -4;
		blast_scores[2][4] = -8;
		blast_scores[3][0] = -4;
		blast_scores[3][1] = -4;
		blast_scores[3][2] = -4;
		blast_scores[3][3] = 5;
		blast_scores[3][4] = -8;
		blast_scores[4][0] = -8;
		blast_scores[4][1] = -8;
		blast_scores[4][2] = -8;
		blast_scores[4][3] = -8;
		blast_scores[4][4] = 0;
	}
	
	private String read_file(String filename) throws IOException {
		InputStream is = new FileInputStream(filename);
		BufferedReader buf = new BufferedReader(new InputStreamReader(is));
		String line = buf.readLine();
		StringBuilder sb = new StringBuilder();
		while (line != null) {
			sb.append(line.trim());
			line = buf.readLine();
		}
		buf.close();
		String fileAsString = sb.toString();
		return fileAsString.replace("\n", "").replace("\r", "");
	}

	public int compute_alignment_score() {
		int score = 0;
		String aligned_s1 = this.sb1.toString();
		String aligned_s2 = this.sb2.toString();
		String aligned_s3 = this.sb3.toString();
		
		for (int i = 0; i < aligned_s1.length(); i++) {
			score += blast_scores[alphabet.indexOf(aligned_s1.charAt(i))][alphabet.indexOf(aligned_s2.charAt(i))] 
					+ blast_scores[alphabet.indexOf(aligned_s1.charAt(i))][alphabet.indexOf(aligned_s3.charAt(i))]
					+ blast_scores[alphabet.indexOf(aligned_s2.charAt(i))][alphabet.indexOf(aligned_s3.charAt(i))];
		}
		return score;
	}

	public int count_matches() {
		int matches = 0;
		String aligned_s1 = this.sb1.toString();
		String aligned_s2 = this.sb2.toString();
		String aligned_s3 = this.sb3.toString();
		
		for (int i = 0; i < aligned_s1.length(); i++) {
			if (aligned_s1.charAt(i) != '-' && aligned_s1.charAt(i) == aligned_s2.charAt(i) && aligned_s1.charAt(i) == aligned_s3.charAt(i))
				matches++;
		}
		return matches;
	}
	
	public int get_aligned_scores()
	{
		int score = 0;
		for (int i = 0; i < this.alignment_scores.size(); i++) {
			score += this.alignment_scores.get(i);
		}
		return score;
	}
	public String get_first_aligned_sequence()
	{
		return this.sb1.toString();
	}
	public String get_second_aligned_sequence()
	{
		return this.sb2.toString();
	}
	public String get_third_aligned_sequence()
	{
		return this.sb3.toString();
	}
	
	/*******************************/
	public static void main(String[] args) throws IOException {
		String currentDir = System.getProperty("user.dir");
		
//		String fn1 = currentDir + "/data/NM_000558.txt";
//		String fn2 = currentDir + "/data/NM_008218.txt";
//		String fn3 = currentDir + "/data/NM_013096.txt";
		
//		String fn1 = currentDir + "/data/NM_010019.txt";
//		String fn2 = currentDir + "/data/NM_001243563.txt";
//		String fn3 = currentDir + "/data/NM_014326.txt";
		
		String fn1 = currentDir + "/data/NM_000545.txt";
		String fn2 = currentDir + "/data/NM_008261.txt";
		String fn3 = currentDir + "/data/NM_000457.txt";
		ThreeDNASequencenAlignment alignment = new ThreeDNASequencenAlignment(fn1, fn2, fn3);
		
		alignment.align();
		int computed_score = alignment.compute_alignment_score();
		int matches = alignment.count_matches();
		
		System.out.println("Computed score: " + computed_score);
		System.out.println("Alignment length: " + alignment.get_first_aligned_sequence().length());
		System.out.println("Exact matches: " + matches);
		System.out.println("Alignment score: " + alignment.get_aligned_scores());
		System.out.println("************");
		System.out.println(alignment.get_first_aligned_sequence());
		System.out.println(alignment.get_second_aligned_sequence());
		System.out.println(alignment.get_third_aligned_sequence());
		System.out.println("************");
		System.out.println("Used Memory: " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) + " bytes");
	}
	/*******************************/
}
