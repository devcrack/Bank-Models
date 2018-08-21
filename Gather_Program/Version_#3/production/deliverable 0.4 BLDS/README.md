# Lanimfe Bank Models
This is a simple suit that allows you get some interesting structure factors like Hards Spheere between others.

## Sinopsis
         ./blds   [--hs]  [--ss] [--yuk] [--exe] [--print] [--plot]
                  [--vf <value>] [--ti <value>]  [--help]
                  
##Descripcion
        Banco de Modelos Lanimfe es un pequeño conjunto de programas fisico/matematicos que permiten realizar calculos para calcular los factores de estructura
        para Esferas Duras, Esfera Suave y Yukawa 3D.

##Opciones
        --hs
            Indica que se esta haciendo referencia al modelo ESFERA DURA


        --ss
            Indica que se esta haciendo referencia al modelo ESFERA SUAVE


        --yuk
            Indica que se esta haciendo referencia al modelo YUKAWA


        --exe 
             Especifica que un calculo se tienen que realizar al modelo que se esta referenciando.
             Se puede referenciar un modelo y realizar un calculo del factor de estructura de dicho modelo. 
            
             El siguiente comando muestra como ejecutaria el calculo para esfera dura. 
           
            ./blds --hs --exe
 
             El siguiente comando muestra como ejecutaria el calculo para esfera dura. 

            ./blds --ss --exe           

             El siguiente comando muestra como ejecutaria el calculo para Yukawa. 

            ./blds --yuk --exe
           
            El siguiente comando muestra como ejecutaria el calculo para Esfera Dura, Esfera Suave y Yukawa, de manera Simultanea. 

            ./blds --hs --hs --yuk --exe


        --vf <value>
            Establece un valor para la Fraccion de Volumen y aplicar dicho valor en el calculo del modelo al que se esta referenciando.
           
            El siguiente comando muestra como realizar la ejecucion del calculo de esfera dura estableciendo un valor para el Factor de Estructura.
           
            ./blds --hs --vf 0.5 --exe
            
            Esto mismo se aplica al resto de los modelos.


        --ti <value>
             Establece un valor para la Temperatura Inicial y aplicar dicho valor en el calculo del modelo al que se esta referenciando.
           
             El siguiente comando muestra como realizar la ejecucion del calculo de esfera dura estableciendo un valor para la Temperatura Incial.
           
            ./blds --hs --ti 0.1 --exe

        
        --plot
            Realiza la ejecucion para la graficacion de los datos del modelo al que se esta referenciando.

            El siguiente comando muestra como realizar la ejecucion para la graficacion de los datos del modelo de Esfera Dura.
           
            ./blds --hs --plot
                  

        --print
           Realiza la ejecucion para la impresion de los datos del modelo al que se esta referenciando.

           El siguiente comando muestra como realizar la ejecucion para la impresion de los datos del modelo de Esfera Dura.           
            ```
           ./blds --hs --plot
           ```

Si no se establece ninguna opcion en la ejecucion del banco de modelos se desplegara el menu interactivo para realizar todas las operaciones previamente mencionadas.
   ```
   ./blds
   ``` 











