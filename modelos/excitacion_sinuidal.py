import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from openseespy.postprocessing.Get_Rendering import plot_model


class ExcitacionSinuidal():

    FILE_DESPLAZAMIENTO_NODO2 = "DESPLAZAMIENTO_NODO2.out"          # DESPLAZAMIENTO N2
    FILE_DESPLAZAMIENTO_NODO1 = "DESPLAZAMIENTO_NODO1.out"          # DESPLAZAMIENTO N1
    FILE_REAC_BASE = "REAC_BASE.out"                                # REACCION_BASE      
    FILE_FUERZAS_COL_SEC = "Fuerzas_Col_sec.out"                    # FUERZAS COL SEC
    FILE_FUERZAS_COL = "FUERZAS_COL.out"                            # FUERZAS COL

    FILE_DEFORMACION_COL = "DEFORMACION_COL.out"                    # DEFORMACION_COL
    FILE_FORCE_COLUMNA = "Data/ForceColSecnumIntgrPts.out"          # FORCE_COLUMNA      
    FILE_DEFORMACION_COLUMNA = "Data/DefoColSec$numIntgrPts.out"    # DEFORMACION_COLUMNA
    FILE_ROTACION_PLASTICA = "PlasticRotation.out"                  # ROTACION PLASTICA     
    FILE_ELEMENTO_AB = "ElementoAB.out"                             # ELEMENTO AB
    GRAVEDAD=9.81


    def __init__(self,nombre='BasicBuilder',numero_dimensiones=2, numero_dots=3, altura_columna=4, fuerza_aplicada=15000 , resistencia_fc=21000 ,dim_travs_colum_altura=1.5, dim_trav_colum_base=1.5, calculo='*'):
        
        self.numero_dimensiones = numero_dimensiones
        self.numero_dots = numero_dots
        self.nombre = nombre
        self.dim_travs_colum_altura = dim_travs_colum_altura
        self.dim_trav_colum_base = dim_trav_colum_base
        self.calculo = calculo
        self.altura_columna = altura_columna
        self.fuerza_aplicada = fuerza_aplicada
        self.resistencia_fc = resistencia_fc
        

        ops.wipe()
        #Creación del modelo 
        ops.model(self.nombre, "-ndm", self.numero_dimensiones, "-ndf", self.numero_dots)

    
    def construir_modelo(self):
        
        #Datos de la estructura
        H_col=1.5 #float(input("Ingrese las dimensiones tranversal de la columna, altura: "))
        B_col=1.5 #float(input("Ingrese las dimensiones tranversal de la columna, base: "))
        h=12 #float(input("Ingrese la altura de la columna: "))
        fuer_apli=8000 #float(input("Ingrese la fuerza aplicada a la columna: "))
        fc=21000 #float(input("Ingrese la resistencia a la compresión: "))
        pes_esp_hor=23.56
        peso_ele=pes_esp_hor*H_col*B_col*h
        peso_estru=peso_ele+fuer_apli
        g=9.81
        masa=peso_estru/g
        masa_1=masa-323.896
        Area=H_col*B_col
        I=(1/12)*B_col*H_col**3
        E=20000000

        #Definir los nodos de la estructura

        ops.node(1, 0, 0)
        ops.node(2, 0, h)

        #Restricciones    
        ops.fix(1, 1, 1, 1)


        #Definir la masa en la cabeza de la columna
        ops.mass(2, masa, 1e-8, 1e-8)

        # resistencia nominal a la compresión del hormigón
        fc1=fc
        Ec=E

        # Definir la sección de columna
        My_col=130000
        Phi_Y_col=0.65e-4
        EI_col_crack=My_col/Phi_Y_col
        b=0.01

        ops.uniaxialMaterial("Steel01", 2, My_col, EI_col_crack, b)
        ops.uniaxialMaterial("Elastic", 3, Ec*Area)
        ops.section("Aggregator", 1, 3, "P", 2, "Mz")

        #Definir la transformacion de coordenadas
        ops.geomTransf("PDelta", 1)
        num_intgr_pts=5
        ops.element("nonlinearBeamColumn", 1, 1, 2, num_intgr_pts, 1, 1)

        # Propiedades del hormigón y acero
        fc1=fc
        Ec=E
        fy=420000
        fyh=420000
        Es=200000000
        shr=0.003

        #Columnas
        dia_bar_long=16
        diametro=dia_bar_long/1000
        num_arr_abajo=8
        num_izq_der=4
        dia_est=10
        estribo=dia_est/1000
        num_ram_Y=2
        num_ram_X=2
        recubrimiento=4
        rec=recubrimiento/100
        sep_est=0.1

        # Definir el modelo de mander
        eco=0.002
        epss=0.004
        esm=0.1
        sp=sep_est*estribo
        rs=(num_ram_X*estribo**2*(np.pi/4)*(B_col-2*rec)+num_ram_Y*estribo**2*(np.pi/4)*(H_col-2*rec))/((B_col-2*rec)*(H_col-2*rec)*sep_est)
        area_var=diametro**2*(np.pi/4)
        area_con=(B_col-2*rec)*(H_col-2*rec)
        rcc=area_var*(num_arr_abajo+num_izq_der)/area_con
        wi=0.2
        ke=((1-wi**2/(6*area_con))*(1-sep_est/(2*B_col))*(1-sep_est/(2*H_col)))/(1-rcc)
        fpl=1/(2*ke*rs*fyh)
        fpcc=(-1.254+2.254*(1+7.94*fpl/fc1)**0.5-2*fpl/fc1)*fc1
        ecc=eco*(1+5*(fpcc/(fc1-1)))
        Esec=fpcc/ecc
        r=Ec/(Ec-Esec)
        ecu=1.5*(0.004+1.4*rs*fyh*(esm/fpcc))
        x=ecu/ecc
        fcu=fpcc*x*r/(r-1+x**r)

        #Definir los materiales de la columna
        #Concreto incofinado
        ops.uniaxialMaterial("Concrete01", 4, -fc1, -eco, 0, -epss)

        #Acero longitudinal
        ops.uniaxialMaterial("Steel01", 5, fy, Es, shr)

        #Concreto cofinado
        ops.uniaxialMaterial("Concrete01", 6, -fpcc, -ecc, -fcu, -ecu)
        ops.uniaxialMaterial("Elastic", 7, Ec)

        #Creacion de las secciones
        y0=B_col/2
        y1=y0-rec
        z0=H_col/2
        z1=z0-rec
        num_div_Z=16
        num_div_Y=16
        num_div_rec=2
        ops.section("Fiber", 10)

        # Definir los parches para recubriento lateral

        crdsI=[-y0, -z0]
        crdsJ=[-y1, -z0]
        crdsK=[-y1, +z0]
        crdsL=[-y0, +z0]

        ops.patch("quad", 4, num_div_rec, (num_div_Z+2*num_div_rec), *crdsI, *crdsJ, *crdsK, *crdsL)

        crdsI=[y1, -z0]
        crdsJ=[y0, -z0]
        crdsK=[y0, +z0]
        crdsL=[y1, +z0]

        ops.patch("quad", 4, num_div_rec, (num_div_Z+2*num_div_rec), *crdsI, *crdsJ, *crdsK, *crdsL)

        # Definir los parches para recubriento inferior y superior
        crdsI=[-y1, -z0]
        crdsJ=[y1, -z0]
        crdsK=[y1, -z1]
        crdsL=[-y1, -z1]

        ops.patch("quad", 4, num_div_Y, num_div_rec, *crdsI, *crdsJ, *crdsK, *crdsL)

        crdsI=[-y1, -z1]
        crdsJ=[y1, -z1]
        crdsK=[y1, +z0]
        crdsL=[-y1, +z0]

        ops.patch("quad", 4, num_div_Y, num_div_rec, *crdsI, *crdsJ, *crdsK, *crdsL)

        # Definir los parches para hormigon confinado
        crdsI=[-y1, -z1]
        crdsJ=[y1, -z1]
        crdsK=[y1, z1]
        crdsL=[-y1, z1]

        ops.patch("quad", 6, num_div_Y, num_div_rec, *crdsI, *crdsJ, *crdsK, *crdsL)

        # Acero superior e inferior
        ops.fiber(0, z1, ((num_arr_abajo/2)*area_var), 5)
        ops.fiber(0, -z1, ((num_arr_abajo/2)*area_var), 5)

        # Acero superior e inferior
        ops.fiber(y1, 0, ((num_izq_der/2)*area_var), 5)
        ops.fiber(y1, 0, ((num_izq_der/2)*area_var), 5)

        ops.section("Aggregator", 2, 7, "Vy", "-section", 10)

        #Definir los recorders
        ops.recorder("Node", "-file", self.FILE_DESPLAZAMIENTO_NODO2, "-time", "-node", 2, "-dof", 1, 2, 3, "disp")
        ops.recorder("Node", "-file", self.FILE_DESPLAZAMIENTO_NODO1, "-time", "-node", 1, "-dof", 1, 2, 3, "disp")
        ops.recorder("Node", "-file", self.FILE_REAC_BASE, "-time", "-node", 1, "-dof", 1, 2, 3, "reaction")

        ops.recorder("Element", "-file", self.FILE_FUERZAS_COL, "-time", "-ele", 1, "globalForce")        
        ops.recorder("Element", "-file", self.FILE_FUERZAS_COL_SEC, "-time", "-ele", 1, "section", 10, "force")

        ops.recorder("Element", "-file", "Defor_col_sec.out", "-time", "-ele", 1, "section", 10, "deformation")
        ops.recorder("Element", "-file", "Fuerzas_Col_sec_num_intgr_pts.out", "-time", "-ele", 1, "section", num_intgr_pts, "force")
        ops.recorder("Element", "-file", "Defor_Col_sec_num_intgr_pts.out", "-time", "-ele", 1, "section", num_intgr_pts, "force")

        #Definir las actuantes sobre la columna
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", -masa)

        #Parametros de analisis estatico
        ops.constraints("Transformation")
        ops.numberer("RCM")
        ops.system("BandGeneral")
        ops.test("EnergyIncr", 1e-8, 5, 0)
        ops.algorithm("ModifiedNewton")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")

        #Analisis modal
        Num_modos=1 #int(input("Ingrese el numero de modos a calcular: "))
        OmegaSq=ops.eigen("-fullGenLapack", Num_modos)
        OmegaSq=np.array(OmegaSq)
        Omega=OmegaSq**0.5
        Periodo=(2*np.pi)/Omega

        for i in range(Num_modos):
            print("Modo", i+1, "T= ", Periodo[i], "s")
            
        # Definir la excitación de la onda sinusoidal
        Tmax=100
        dt=0.01
        g=9.81
        Ampl_ond_sinusoidal=0.5*g
        T_onda_sinusoidal=0.35
        Dur_onda_sinusoidal=5
        Dir=1
        Npuntos=1000

        Omega_onda_sinusoidal=(2*np.pi/T_onda_sinusoidal)
        Vel0=Ampl_ond_sinusoidal*(-1)/Omega_onda_sinusoidal
        ops.timeSeries("Sine", 0, Dur_onda_sinusoidal, T_onda_sinusoidal, "-factor", Ampl_ond_sinusoidal)
        ops.pattern("UniformExcitation", 8, Dir, "-accel", 0, "-vel0", Vel0)

        #Parametros de analisis estatico

        ops.constraints("Plain")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("NormDispIncr", 1.0e-8, 10)
        ops.algorithm("Newton")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        ops.analyze(Npuntos, dt)

        #ops.remove("recorders")

        if self.calculo == 'REACCION_BASE' :                
            array_mek=np.loadtxt(self.FILE_REAC_BASE)
            self.construye_figura(array_mek,serie=[1,2,3], recoder=self.calculo)                            
        elif self.calculo == 'DESPLAZAMIENTO_NODO2':            
            array_mek=np.loadtxt(self.FILE_DESPLAZAMIENTO_NODO2)
            self.construye_figura(array_mek,serie=[1,2,3], recoder=self.calculo)
        elif self.calculo == 'DESPLAZAMIENTO_NODO1':            
            array_mek=np.loadtxt(self.FILE_DESPLAZAMIENTO_NODO1)
            self.construye_figura(array_mek,serie=[1,2,3], recoder=self.calculo)
        
        elif self.calculo == 'FUERZAS_COL':            
            array_mek=np.loadtxt(self.FILE_FUERZAS_COL)           
            self.construye_figura(array_mek,serie=[1], recoder=self.calculo)

        elif self.calculo == 'FUERZAS_COL_N1':            
            array_mek=np.loadtxt(self.FILE_FUERZAS_COL_N1)
            self.construye_figura(array_mek,serie=[1], recoder=self.calculo)
        
        elif self.calculo == 'DEFORMACION_COL':            
            array_mek=np.loadtxt(self.FILE_DEFORMACION_COL)
            self.construye_figura(array_mek,serie=[1], recoder=self.calculo)

        elif self.calculo == 'FORCE_COLUMNA':            
            array_mek=np.loadtxt(self.FILE_FORCE_COLUMNA)
            self.construye_figura(array_mek,serie=[1], recoder=self.calculo)

        elif self.calculo == 'DEFORMACION_COLUMNA':            
            array_mek=np.loadtxt(self.FILE_DEFORMACION_COLUMNA)
            self.construye_figura(array_mek,serie=[1], recoder=self.calculo)

        else:            
            
            #datos_1=np.loadtxt(self.FILE_REAC_BASE)
            #self.construye_figura(array_mek,serie=[1,2,3], recoder='REACCION_BASE')

            datos_2=np.loadtxt(self.FILE_DESPLAZAMIENTO_NODO2)
            #self.construye_figura(array_mek,serie=[1,2,3], recoder='DESPLAZAMIENTO_NODO2')

            #datos_3=np.loadtxt(self.FILE_DESPLAZAMIENTO_NODO1)
            #self.construye_figura(array_mek,serie=[1,2,3], recoder='DESPLAZAMIENTO_NODO1')

            """array_mek=np.loadtxt(self.FILE_FUERZAS_COL)
            self.construye_figura(array_mek,serie=[1], recoder='FUERZAS_COL')

            array_mek=np.loadtxt(self.FILE_FUERZAS_COL_N1)
            self.construye_figura(array_mek,serie=[1], recoder='FUERZAS_COL_N1')

            array_mek=np.loadtxt(self.FILE_DEFORMACION_COL)
            self.construye_figura(array_mek,serie=[1], recoder='DEFORMACION_COL')

            array_mek=np.loadtxt(self.FILE_FORCE_COLUMNA)
            self.construye_figura(array_mek,serie=[1], recoder='FORCE_COLUMNA')

            array_mek=np.loadtxt(self.FILE_DEFORMACION_COLUMNA)
            self.construye_figura(array_mek,serie=[1], recoder='DEFORMACION_COLUMNA')"""

            serie_total = [ [ datos_2, [1,2,3],'DESPLAZAMIENTO_NODO2' ] ]

            if serie_total:
                self.construye_figura(serie_total=serie_total, recoder='*')
            

    def concrete_section(self, core_tag, cover_tag, B, H, cover, nfcore, nfcover):
        ops.patch("rect", core_tag, nfcore, 1, -H/2+cover, -B/2+cover, H/2+cover, B/2+cover)
        ops.patch("rect", cover_tag, nfcover, 1, -H/2, -B/2, -H/2+cover, B/2)
        ops.patch("rect", cover_tag, nfcover, 1, H/2-cover, -B/2, H/2, B/2)
        ops.patch("rect", cover_tag, nfcore, 1, -H/2+cover, -B/2, H/2-cover, B/2+cover)
        ops.patch("rect", cover_tag, nfcore, 1, -H/2+cover, B/2-cover, H/2-cover, B/2)

    
    def construye_figura(self, array_mek=None, serie=[], recoder=None, serie_total=[], resultados=[]):
    
        """Construye una figura basado en un array de datos, una serie de gráficos        
        """
        
        if recoder != '*':

            nombre_figura = "Gráfico # "
            titulo = None

            for i in serie:

                if recoder != None:
                    titulo = recoder+' : '+nombre_figura+str(i)
                    plt.figure(titulo)
                else:
                    titulo = nombre_figura + str(i)
                    plt.figure(titulo)
                
                plt.plot(array_mek[:,0], array_mek[:,i])
                
                plt.ylabel('Desplazamiento horizontal (m)')
                plt.xlabel('Tiempo (s)')
                plt.show()
        else :
            
            from PIL import Image
            import matplotlib.image as mpimg
            import os
            from pathlib import Path
            
            BASE_DIR = Path(__file__).resolve().parent.parent

            IMAGEN_COLUMNA = os.path.join(os.path.join(BASE_DIR,'modelos'), 'columna.jpeg')

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
            ax1.imshow(pil_img) 
            ax1.axes.xaxis.set_visible(False)
            ax1.axes.yaxis.set_visible(False)

            ax2 = plt.subplot2grid((FILAS, COLUMNAS), (2, 0))
            ax2.axes.xaxis.set_visible(False)
            ax2.axes.yaxis.set_visible(False)
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)
            ax2.spines['bottom'].set_visible(False)
            ax2.spines['left'].set_visible(False)  

            ax6 = plt.subplot2grid((FILAS, COLUMNAS), (0, 1))
            ax7 = plt.subplot2grid((FILAS, COLUMNAS), (1, 1))
            ax8 = plt.subplot2grid((FILAS, COLUMNAS), (2, 1))

                                                
            for columna, elemento in enumerate(serie_total):                
                
                array_mek = elemento[0]
                serie = elemento[1]
                calculo_var = elemento[2]

                if calculo_var == 'DESPLAZAMIENTO_NODO2':
                                        
                    ax6.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:red')
                    ax6.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=8)
                    ax6.set_xlabel('Tiempo (s)')
                    ax6.set_ylabel('Desplazamiento H. (m)')
                    
                    ax7.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:red')
                    ax7.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}", fontsize=8)
                    ax7.set_xlabel('Tiempo (s)')
                    ax7.set_ylabel('Desplazamiento H. (m)')

                    ax8.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:red')
                    ax8.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=8)
                    ax8.set_xlabel('Tiempo (s)')
                    ax8.set_ylabel('Desplazamiento H. (m)')      


                        
            fig.suptitle("MODELOS MATEMÁTICOS: EXCITACION SINUIDAL")                             
            fig.subplots_adjust(left=0.305, bottom=0.076, right=0.92, top=0.88, wspace=0.430, hspace=1)
    
        mng = plt.get_current_fig_manager()
        mng.set_window_title("MODELOS MATEMÁTICOS UTPL")
        mng.resize(*mng.window.maxsize())                
        plt.show() 

    def mostrar_resultados(self, lista_mensajes=[],ax=None):
                
        fig, ax1 = plt.subplots(1, 1, constrained_layout=True, sharey=True)       
        ax1.set_xlabel('Tiempo (s)')
        ax1.set_ylabel('Desplazamiento H. (m)')   
                        
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

    
        