# HE.CNNtraining

## Privacy-Preserving CNN Training with Transfer Learning


### Sigmoid(x) ~ poly3 = 0.5 + 0.10679534503216294.*x  + -0.00038503259805075.*x.^3;  (lambda = 128)


HE.CNNtraining is a project for implementing our CNN training on  encrypted MNIST images (Privacy-Preserving CNN Training with Transfer Learning)

## How to run this program? 

### Dependencies

In a Ubuntu cloud, our implementation requires the following libraries in order:
* `g++`:      
```sh
               # apt install g++ 
```

* `make`:       
```sh
                # apt install make
```

* `m4`: #        
```sh
                 # apt install m4
```

* `GMP`(ver. 6.1.2):      
```sh
                           # cd gmp-x.x.x  
                           # ./configure --enable-cxx  
                           # make
                           # make install
                           # ldconfig
```

* `NTL`(ver. 11.3.0): 
```sh
                     # cd ntl-x.x.x
                     # cd src
                     # ./configure NTL_THREADS=on NTL_THREAD_BOOST=on NTL_EXCEPTIONS=on
                     # make
                     # make install
```

### Running CNNinference 

You need to configure and build the HE.CNNtraining project. If on a Ubuntu 22.04 x64 you placed the project in the path:
```sh
/home/john/eclipse-workspace/MyHECNNtraining/$ls
data  Debug  Default  HEAAN  lib  result  run  src
```
It might be much easy to configure and build the project.  

After that, in the 'Default' folder, you can run our project by the following command lines:

```sh
# make clean
# make all
# ./MyHECNNtraining
``` 

You can change the source codes and then repeat the above lines to debug your own project.

## Running a test source code

In the 'Default' folder, you can find one running results:   

        '2023Apr07with12vCPUs_nohup.out'  
        

     
         


            
            
    

