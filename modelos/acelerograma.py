import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from openseespy.postprocessing.Get_Rendering import plot_model


class Acelerograma():

    #Variables de Clase
    area = None
    numero_pasos=10
    E=20000000
    G=9.81

    #RECODERS
    FILE_DESPLAZAMIENTO = "DESPLAZAMIENTO.out"  # DESPLAZAMIENTO
    FILE_REAC_BASE = "REAC_BASE.out"            # REACCION_BASE
    FILE_SOLIC_COL = "SOLIC_COL.out"            # SOLICITACION

    def __init__(self,nombre='BasicBuilder',numero_dimensiones=2, numero_dots=3, data_file='datos.txt', altura_columna=12, masa=1500, dim_travs_colum_altura=1.5, dim_trav_colum_base=1.5, calculo='*'):

        self.nombre = nombre
        self.numero_dimensiones = numero_dimensiones
        self.numero_dots = numero_dots
        self.data_file = data_file

        #Datos de la Estructura
        self.altura_columna = altura_columna
        self.masa = masa
        self.dim_travs_colum_altura = dim_travs_colum_altura
        self.dim_trav_colum_base = dim_trav_colum_base
        self.calculo = calculo

        ops.wipe()
        ops.model(self.nombre, "-ndm", self.numero_dimensiones, "-ndf", self.numero_dots)

    def construir_modelo(self):

        """Construye el modelo con los valores proporcionados por el usuario en el constructor de clase
        """

        area=self.dim_travs_colum_altura*self.dim_trav_colum_base
        I=(1/12)*self.dim_trav_colum_base*self.dim_travs_colum_altura**3

        #Definir los nodos de la estructura
        ops.node(1, 0, 0)
        ops.node(2, 0, self.altura_columna)

        #Definir la masa de la estructura
        ops.mass(2, self.masa, 0, 0)

        #Definir las restricciones
        ops.fix(1, 1, 1, 1)

        #Definir la transformacion geometrica
        ops.geomTransf("Linear", 1)

        #Definir elemento de la columna
        ops.element("elasticBeamColumn", 1, 1, 2, area, self.E, I, 1)

        #Definir los recorders
        ops.recorder("Node", "-file", self.FILE_DESPLAZAMIENTO, "-time", "-node", 2, "-dof", 1, 2, 3, "disp")
        ops.recorder("Node", "-file", self.FILE_REAC_BASE, "-time", "-node", 1, "-dof", 1, 2, 3, "reaction")
        ops.recorder("Element", "-file", self.FILE_SOLIC_COL, "-time", "-ele", 1, "force", 1)

        #Definir las actuantes sobre la columna
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 0, (-1*self.masa*self.G), 0)

        #Parametros de analisis estatico

        ops.constraints("Plain")
        ops.numberer("RCM")
        ops.system("BandGeneral")
        ops.test("NormDispIncr", 1.0e-6, 10)
        ops.algorithm("Newton")        
        ops.integrator("LoadControl", 1/self.numero_pasos)
        ops.analysis("Static")
        ops.analyze(10)

        #Definir análisis dinamico
        ops.loadConst("-time", 0)
        numero_nodos=1
        factor=0.01
        dt=0.01
        numero_puntos=10000

        valor_dir=1
        ops.timeSeries("Path", 2, "-dt", dt, "-filePath", self.data_file, "-factor", factor)
        ops.pattern("UniformExcitation", 2, valor_dir, "-accel", 2)

        #Parametros de analisis estatico

        ops.constraints("Plain")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("NormDispIncr", 1.0e-8, 10)
        ops.algorithm("Newton")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        ops.analyze(numero_puntos, dt)

        #Definir el amortiguamiento
        omegasq=ops.eigen("-fullGenLapack", 1)
        omegasq=np.array(omegasq)
        omega=omegasq**0.5
        periodo=(2*np.pi)/omega
        ops.rayleigh(0, 0, 0, 0.01)

        for i in range(numero_nodos):
            print("Modo", i+1, "T= ", periodo[i], "s")
            
        ops.remove("recorders")
            
        #Personalizar la carga de los datos
        # 
        serie_total= []        

        if self.calculo == 'REACCION_BASE' :                
            array_mek=np.loadtxt(self.FILE_REAC_BASE)
            self.construye_figura(array_mek,serie=[1,2,3], recoder=self.calculo)                            

        elif self.calculo == 'SOLICITACION':            
            array_mek=np.loadtxt(self.FILE_SOLIC_COL)
            self.construye_figura(array_mek,serie=[1], recoder=self.calculo)

        elif self.calculo == 'DESPLAZAMIENTO':            
            array_mek=np.loadtxt(self.FILE_DESPLAZAMIENTO)
            self.construye_figura(array_mek,serie=[1,2,3], recoder=self.calculo)

        else:            
            datos_1=np.loadtxt(self.FILE_REAC_BASE)                        
            datos_2=np.loadtxt(self.FILE_SOLIC_COL)                        
            datos_3=np.loadtxt(self.FILE_DESPLAZAMIENTO)

            serie_total = [[ datos_1, [1,2,3],'REACCION_BASE' ], [ datos_2, [1],    'SOLICITACION' ], [ datos_3, [1,2,3],'DESPLAZAMIENTO' ] ]
            #self.construye_figura(serie_total=serie_total, recoder='*')
        
        nodo=[1, 2]
        dof=[1, 2, 3]

        nodo1=nodo[0]
        nodo2=nodo[1]

        gra_liber1=dof[0]
        gra_liber2=dof[1]
        gra_liber3=dof[2]

        #Resultados

        RESULTADOS = []

        RESULTADOS.append(str("El desplazamiento en el nodo: " + str(nodo2)+ " para el grado de libertad: "+str(gra_liber1)+" es: "+str(ops.nodeDisp(2, 1))+" m") )        
        RESULTADOS.append(str("El desplazamiento en el nodo: " + str(nodo2)+ " para el grado de libertad: "+str(gra_liber2)+" es: "+str(ops.nodeDisp(2, 2))+" m") )        
        RESULTADOS.append(str("El desplazamiento en el nodo: " + str(nodo2)+ " para el grado de libertad: "+str(gra_liber3)+" es: "+str(ops.nodeDisp(2, 3))+" m") )        

        RESULTADOS.append(str("La reaccción en el nodo: " + str(nodo1)+ " para el grado de libertad: "+str(gra_liber1)+" es: "+str(ops.nodeReaction(1, 1))+" N" ))
        RESULTADOS.append(str("La reaccción en el nodo: " + str(nodo1)+ " para el grado de libertad: "+str(gra_liber2)+" es: "+str(ops.nodeReaction(1, 2))+" N") )
        RESULTADOS.append(str("La reaccción en el nodo: " + str(nodo1)+ " para el grado de libertad: "+str(gra_liber3)+" es: "+str(ops.nodeReaction(1, 3))+" N") )

        RESULTADOS.append(str("La reaccción en el elemento: "+ str(nodo1)+ " para el grado de libertad:"+ str(gra_liber1)+" es:"+str(ops.nodeReaction(1, 1))+" N"))

        if serie_total:
            self.construye_figura(serie_total=serie_total, recoder='*', resultados=RESULTADOS)

        self.mostrar_resultados(lista_mensajes=RESULTADOS)

        #plot_model("nodes","elements") 




    def construye_figura(self, array_mek=None, serie=[], recoder=None, serie_total=[], resultados=[]):
        
        """Construye una figura basado en un array de datos, una serie de gráficos        
        """
                
        titulo = None        
        
        if recoder != '*':

            ventana_principal = plt.figure(0)
            ventana_principal.canvas.set_window_title(recoder)

            for i in serie:
                plt.subplots_adjust(hspace=0.70)
                plt.subplot(2,2,i)
                titulo = "RESPUESTA : GRADO DE LIBERTAD #{}".format(i)            
                plt.title(titulo)                        
                plt.xlabel('Tiempo (seg)')
                plt.ylabel('Desplazamiento (m)')
                print ("Normal :",i)
                plt.plot(array_mek[:,0], array_mek[:,i],label="Resultado")
                plt.legend()                           
        else :
            
            from PIL import Image
            import matplotlib.image as mpimg
            import os
            from pathlib import Path
            
            BASE_DIR = Path(__file__).resolve().parent.parent

            IMAGEN_COLUMNA = os.path.join(os.path.join(BASE_DIR,'modelos'), 'Columna_con_carga_ciclica.jpeg')

            pil_img = Image.open(IMAGEN_COLUMNA)
            

            num_subplots = 0
            
            for elemento in serie_total:
                num_subplots += len(elemento[1]) 
            
            if num_subplots % 2 != 0:
                num_subplots+=1
                        
            fig = plt.figure()


            FILAS = 3
            COLUMNAS = 3
            ax1 = plt.subplot2grid((FILAS, COLUMNAS), (0, 0), rowspan=2)
            ax2 = plt.subplot2grid((FILAS, COLUMNAS), (2, 0))
            ax3 = plt.subplot2grid((FILAS, COLUMNAS), (0, 1))
            ax4 = plt.subplot2grid((FILAS, COLUMNAS), (1, 1))
            ax5 = plt.subplot2grid((FILAS, COLUMNAS), (2, 1))
            ax6 = plt.subplot2grid((FILAS, COLUMNAS), (0, 2))
            ax7 = plt.subplot2grid((FILAS, COLUMNAS), (1, 2))
            ax8 = plt.subplot2grid((FILAS, COLUMNAS), (2, 2))
            
            #Para controlar el subplot de Resultados se visualice
            if COLUMNAS == 4: 
                ax9 = plt.subplot2grid((3,4), (0, 3), rowspan=3)
            

            for columna, elemento in enumerate(serie_total):                
                
                array_mek = elemento[0]
                serie = elemento[1]
                calculo_var = elemento[2]
            
                if calculo_var == 'REACCION_BASE':  
                    
                    ax3.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:orange')
                    ax3.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=8)
                    ax3.set_xlabel('Tiempo (seg)')
                    ax3.set_ylabel('Desplazamiento (m)')
                                        
                    ax4.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:orange')
                    ax4.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}", fontsize=8)
                    ax4.set_xlabel('Tiempo (seg)')
                    ax4.set_ylabel('Desplazamiento (m)')
                    
                    ax5.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:orange')
                    ax5.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=8)
                    ax5.set_xlabel('Tiempo (seg)')
                    ax5.set_ylabel('Desplazamiento (m)')

                elif calculo_var == 'SOLICITACION':
                                        
                    ax1.imshow(pil_img) 
                    ax1.axes.xaxis.set_visible(False)
                    ax1.axes.yaxis.set_visible(False)
                    
                    ax2.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:green')  
                    ax2.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=8)
                    ax2.set_xlabel('Tiempo (seg)')
                    ax2.set_ylabel('Desplazamiento (m)')

                else:

                    ax6.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:red')
                    ax6.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=8)
                    ax6.set_xlabel('Tiempo (seg)')
                    ax6.set_ylabel('Desplazamiento (m)')
                    
                    ax7.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:red')
                    ax7.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}", fontsize=8)
                    ax7.set_xlabel('Tiempo (seg)')
                    ax7.set_ylabel('Desplazamiento (m)')

                    ax8.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:red')
                    ax8.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=8)
                    ax8.set_xlabel('Tiempo (seg)')
                    ax8.set_ylabel('Desplazamiento (m)')                
                        
            if COLUMNAS == 4:
                self.mostrar_resultados(lista_mensajes=resultados,ax=ax9)

            fig.suptitle("MODELOS MATEMÁTICOS: COLUMNA DE HORMIGÓN ARMADO CON MOVIMIENTO DE BASE")                             
            fig.subplots_adjust(left=0.05, bottom=0.076, right=0.971, top=0.88, wspace=0.562, hspace=1)
    
        mng = plt.get_current_fig_manager()
        mng.set_window_title("MODELOS MATEMÁTICOS UTPL")
        mng.resize(*mng.window.maxsize())                
        plt.show() 

    def mostrar_resultados(self, lista_mensajes=[],ax=None):
                
        fig, ax1 = plt.subplots(1, 1, constrained_layout=True, sharey=True)       
        ax1.set_xlabel('Tiempo (seg)')
        ax1.set_ylabel('Desplazamiento (m)')   
                        
        resultado = ''
        
        for index, cadena in enumerate(lista_mensajes):           
            resultado += cadena+'\n'
                        
        ax1.axes.xaxis.set_visible(False)
        ax1.axes.yaxis.set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.text(0.038, 0.75, resultado , horizontalalignment='left', verticalalignment='center_baseline', transform=ax1.transAxes, bbox={'facecolor': 'orange', 'alpha': 0.25, 'pad': 5}, fontsize=9)        
        fig.suptitle('\n\nResultados', fontsize=16)

        plt.show()

        
        

        

        

