import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt


class Momento_Curvatura():

    #RECODERS
    FILE_MOMENTO_CURVATURA = "Momento_curvatura.out"  # MOMENTO CURVATURA
    FILE_MATERIAL_1_TOP = "Material1+Top.out"  # MATERIAL 1 TOP
    FILE_MATERIAL_2_TOP = "Material2+Top.out"  # MATERIAL 2 TOP
    
    FILE_MATERIAL_1_BOT = "Material1+Bot.out"  # MATERIAL 1 BOT
    FILE_MATERIAL_2_BOT = "Material2+Bot.out"  # MATERIAL 2 BOT

    
    def __init__(self,nombre='BasicBuilder',numero_dimensiones=2, numero_dots=3, altura_columna=1.842, dim_travs_colum_altura=0.50, dim_trav_colum_base=0.50, calculo='*'):
                
        self.nombre = nombre
        self.numero_dimensiones = numero_dimensiones
        self.numero_dots = numero_dots                

        #Datos de la Estructura
        self.altura_columna = altura_columna        
        self.dim_travs_colum_altura = dim_travs_colum_altura
        self.dim_trav_colum_base = dim_trav_colum_base
        self.calculo = calculo
        
        ops.wipe()
        ops.model("basic", "-ndm", self.numero_dimensiones, "-ndf", self.numero_dots)

        
    def construir_modelo(self):

        #Datos
        self.dim_travs_colum_altura=0.50 #float(input("Ingrese las dimensiones tranversal de la columna, altura: "))
        self.dim_trav_colum_base=0.50 #loat(input("Ingrese las dimensiones tranversal de la columna, base: "))
        #h=1.842 #t(input("Ingrese la altura de la columna: "))

        fc=28000
        rec=0.04
        dia_varcol=16        
        var_sup_inf=4
        var_der_izq=4
        areaFiber=((np.pi*dia_varcol**2)/4)/1000**2

        # Definir los materiales para la columna
        # ------------------------------------------
        # CONCRET0                   
        # Concreto confinado
        fpc=-fc
        epsc0=-0.003
        fpcu=0.2*fpc
        epsU=-0.006
        lambdaC=0.5
        ft=(0.62*(28)**0.5)*1000
        Ets=(2*(fpc)/epsc0)

        ops.uniaxialMaterial("Concrete02", 1 , fpc , epsc0 , fpcu , epsU , lambdaC , ft , Ets)

        # ACERO DE REFUERZO
        # Acero de refuerzo
        Fy=420000
        E0=200000000
        b=0.007

        ops.uniaxialMaterial( 'Steel02' , 2 , Fy , E0 , b)

        #Definir la sección a la que se va a utilizar
        GJ=1e14
        ops.section( 'Fiber' , 1 , '-GJ' , GJ )

        numSubdivY=20
        numSubdivZ=1
        B=0.40
        H=0.40
        crdsI=[-H/2, -B/2]
        crdsJ=[H/2, B/2]

        ops.patch( 'rect' , 1 , numSubdivY , numSubdivZ , * crdsI , * crdsJ )

        #Definir el acero en la viga
        numFiberSup=var_sup_inf
        areaFiber6=284e-6
        inicio=[-H/2+0.05, -B/2+0.05]
        fin=[-H/2+0.05, B/2-0.05]

        ops.layer( 'straight' , 2 , numFiberSup , areaFiber , * inicio , * fin )

        numFiberInf=var_sup_inf
        areaFiber7=387e-6
        inicio=[H/2-0.05, -B/2+0.05]
        fin=[H/2-0.05, B/2-0.05]

        ops.layer( 'straight' , 2 , numFiberInf , areaFiber , * inicio , * fin )

        numFiberSup=var_der_izq
        areaFiber6=284e-6
        inicio=[H/2+0.05, B/2+0.05]
        fin=[-H/2+0.05, B/2-0.05]

        ops.layer( 'straight' , 2 , numFiberSup , areaFiber , * inicio , * fin )

        numFiberInf=var_der_izq
        areaFiber7=387e-6
        inicio=[H/2-0.05, -B/2+0.05]
        fin=[-H/2-0.05, -B/2-0.05]

        ops.layer( 'straight' , 2 , numFiberInf , areaFiber , * inicio , * fin )

        # Analisis de momento-curvatura
        # Carga axial en la sección
        P= 750;

        # Nodos
        ops.node(1, 0.0, 0.0)
        ops.node(2, 0.0, 0.0)

        # Restriccuiones
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 0, 1, 0)

        # Elemento de longitud cero
        ops.element("zeroLengthSection", 1, 1, 2, 1)

        # Patron de carga axial
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -P, 0.0, 0.0)

        #Crear el analisis
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("BandGeneral")
        ops.test("EnergyIncr", 1.0e-6, 30)
        ops.algorithm("Newton")
        NSteps=10
        ops.integrator("LoadControl", 1/NSteps)
        ops.analysis("Static")

        ok=ops.analyze(NSteps)
        if ok==0:
            print("CARGA AXIAL APLICADA")
        else:
            print("ERROR AL APLICAR LA CARGA AXIAL")
        ops.wipeAnalysis()

        ops.loadConst("-time", 0.0)

        if ok==0:
            # Definir patron de carga
            ops.timeSeries("Linear", 2)
            ops.pattern("Plain", 2, 2)
            ops.load(2, 0.0, 0.0, 1.0)
            
            #Crear el analisis
            ops.constraints("Plain")
            ops.numberer("Plain")
            ops.system("BandGeneral")
            ops.test("EnergyIncr", 1.0e-6, 1000)
            ops.algorithm("Newton")
            maxk=0.0021/(0.7*H)*20
            numIncr=1000
            ops.integrator("DisplacementControl", 2, 3, maxk/numIncr)
            ops.analysis("Static")
            
            # Recorder almacena Momento-Deformacion
            ops.recorder("Node", "-file", self.FILE_MOMENTO_CURVATURA, "-time", "-node", 2,
                        "-dof", 1, 3, "disp")
            
            # Recorder almacena Momento-Deformacion
            ops.recorder("Element", "-file", self.FILE_MATERIAL_1_TOP, "-time","-ele", 
                        1, "section", "fiber", "1.0", "0.0", "1", "stressStrain" )
            
            ops.recorder("Element", "-file", self.FILE_MATERIAL_2_TOP, "-time","-ele", 
                        1, "section", "fiber", "1.0", "0.0", "2", "stressStrain" )
            
            ops.recorder("Element", "-file", self.FILE_MATERIAL_1_BOT, "-time","-ele", 
                        1, "section", "fiber", "-1.0", "0.0", "1", "stressStrain" )
            
            ops.recorder("Element", "-file", self.FILE_MATERIAL_2_BOT, "-time","-ele", 
                        1, "section", "fiber", "-1.0", "0.0", "2", "stressStrain" )
            
            while ops.nodeDisp(2, 3) <= maxk:
                ok=ops.analyze(1)
                if ok != 0:
                    print("EL ANALISIS DE MOMENTO CURVATURA FALLO: ", ops.nodeDisp(2, 3))
                    break
            if ok == 0:
                print("EL ANALISIS TERMINO EXITOSAMENTE")

            ops.wipeAnalysis()
            ops.remove("recorders")
            
        datos_1=np.loadtxt(self.FILE_MOMENTO_CURVATURA)  

        serie_total = [ [ datos_1, [1,2,3],'MOMENTO_CURVATURA'] ]

        if serie_total:
            self.construye_figura(serie_total=serie_total, recoder='*')

    
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
                plt.xlabel('Desplazamiento H. (m)')
                plt.ylabel('Fuerza Cortante (Kn)')
                print ("Normal :",i)
                plt.plot(array_mek[:,0], array_mek[:,i],label="Resultado")
                plt.legend()                           
        else :
            
            from PIL import Image
            import matplotlib.image as mpimg
            import os
            from pathlib import Path
            
            BASE_DIR = Path(__file__).resolve().parent.parent
            IMAGEN_COLUMNA = os.path.join(os.path.join(BASE_DIR,'modelos'), 'columna.jpeg')

            pil_img = Image.open(IMAGEN_COLUMNA)

            IMAGEN_COLUMNA_ENSAYOS = os.path.join(os.path.join(BASE_DIR,'modelos'), 'labEnsayo1.jpg')
            pil_img_ensayos = Image.open(IMAGEN_COLUMNA_ENSAYOS)
            

            num_subplots = 0
            
            for elemento in serie_total:
                num_subplots += len(elemento[1]) 
            
            if num_subplots % 2 != 0:
                num_subplots+=1
                        
            fig = plt.figure()


            FILAS = 2
            COLUMNAS = 3

            ax1 = plt.subplot2grid((FILAS, COLUMNAS), (0, 0), rowspan=2)
            ax1.imshow(pil_img) 
            ax1.axes.xaxis.set_visible(False)
            ax1.axes.yaxis.set_visible(False)

            """ax2 = plt.subplot2grid((FILAS, COLUMNAS), (2, 0))
            ax2.axes.xaxis.set_visible(False)
            ax2.axes.yaxis.set_visible(False)
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)
            ax2.spines['bottom'].set_visible(False)
            ax2.spines['left'].set_visible(False)"""

            ax3 = plt.subplot2grid((FILAS, COLUMNAS), (0, 1), rowspan=2)

            ax6 = plt.subplot2grid((FILAS, COLUMNAS), (0, 2), rowspan=2)
            ax6.imshow(pil_img_ensayos) 
            ax6.axes.xaxis.set_visible(False)
            ax6.axes.yaxis.set_visible(False)

    
            for columna, elemento in enumerate(serie_total):                
                
                array_mek = elemento[0]
                serie = elemento[1]
                calculo_var = elemento[2]
            
                if calculo_var == 'MOMENTO_CURVATURA':
                    ax3.plot(array_mek[:,2], array_mek[:,0],'tab:orange')
                    ax3.set_title(f" {calculo_var}", fontsize=8)
                    ax3.set_xlabel('Desplazamiento H. (m)')
                    ax3.set_ylabel('Fuerza Cortante (Kn)')
                else:
                    pass
            
            fig.suptitle("ENSAYOS MATEMÁTICOS: MOMENTO CURVATURA")                             
            fig.subplots_adjust(left=0.05, bottom=0.076, right=0.971, top=0.88, wspace=0.562, hspace=1)
    
        mng = plt.get_current_fig_manager()
        mng.set_window_title("ENSAYOS MATEMÁTICOS UTPL")
        mng.resize(*mng.window.maxsize())                
        plt.show() 

    def mostrar_resultados(self, lista_mensajes=[],ax=None):
                
        fig, ax1 = plt.subplots(1, 1, constrained_layout=True, sharey=True)       
        ax1.set_xlabel('Desplazamiento H. (m)')
        ax1.set_ylabel('Fuerza Cortante (Kn)')   
                        
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
