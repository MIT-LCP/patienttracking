package com.ikarosilva.pttrack;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;


public class  DistanceMatrix implements Callable<Boolean>{
	private static final int MAX_THREADS=Runtime.getRuntime().availableProcessors();
	private final int N;		//Total number of records
	private final BlockingQueue<Integer> queue;
	private final double[] data;
	private double[][] distance;
	private final int threads;
	private final int step;
	public DistanceMatrix(double[] data, int th) throws InterruptedException{
		//Set task queue 
		N= data.length;
		queue=new ArrayBlockingQueue<Integer>(N);
		this.data=data;
		distance=new double[N][N];
		th=(th < 1) ? MAX_THREADS:th; //In case user enters negative or zero
		threads=(th > MAX_THREADS) ? MAX_THREADS:th;

		//Queue is the starting indices for of data for which each thread will be responsible for
		step=N/threads;
		int index=0;
		for(int n=0;n<threads;n++){
			queue.put(index);
			index+=step;
		}
	}

	public double[][] getDistanceMatrix(){
		return distance;
	}

	public static int getNumberOfProcessors(){return MAX_THREADS;}

	public int getNumberOfThreads(){
		return threads;
	}

	public Boolean call(){
		Integer taskInd;
		//along id=Thread.currentThread().getId();
		while ((taskInd = queue.poll()) != null ){ 
			//System.err.println("Thread [" + id + "]: Processing: " + taskInd);
			compute(taskInd);
		}
		return true;
	}


	private void compute(int start){
		//Computer distances of columns from start to end 
		// for all rows
		int end=((start+step)>(N-1)) ? (N-1):(start+step);
		for(int row=start;row<end;row++){
			for(int col=row+1;col<N;col++){
				distance[row][col]=Math.pow(data[row]-data[col],2);
				distance[col][row]=distance[row][col];
			}
		}
	}

	public void start() throws Exception{
		
		ArrayList<Future<Boolean>> futures=
				new ArrayList<Future<Boolean>>(threads);
		ExecutorService executor= 
				Executors.newFixedThreadPool(threads);
		System.err.println("Starting computation with : " + threads + " threads.");
		try {
			for(int i=0;i<threads;i++){
				Future<Boolean> future= executor.submit(this);
				futures.add(future);
			}
			for( Future<Boolean> future: futures){
				if(future.get()==false)
					System.err.println("Future computation failed:  " + future.toString());
			}
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		} 
		executor.shutdown();
	}


	public static void main(String[] args) throws Exception {
		int N=10;
		double[] data=new double[N];
		for(int n=0;n<N;n++)
			data[n]=Math.random();		
		DistanceMatrix d=new DistanceMatrix(data,0);
		d.start();
		double[][] dist=d.getDistanceMatrix();
		for(int n=0;n<N;n++){
			for(int m=0;m<N;m++){
				System.out.print(" " + dist[n][m] + " ");
			}
			System.out.println("");
		}
	}

}
