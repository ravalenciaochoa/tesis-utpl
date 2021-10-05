import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from openseespy.postprocessing.Get_Rendering import plot_model, plot_modeshape

class CargaDistribuida():

    CARGA_DISTRIBUIDA = 150
    PES_ESP_HOR = 23.56
    GRAVEDAD = 9.81
    
    FILE_DESPLAZAMIENTO = "DESPLAZAMIENTO.out"          # DESPLAZAMIENTO N2
    FILE_REAC_BASE = "REAC_BASE.out"                    # REACCION_BASE 


    def __init__(self, nombre='BasicBuilder', numero_dimensiones=2, numero_dots=3, altura_columna=12, fuerza_aplicada=500, resistencia_comprension=21000 ,dim_travs_colum_altura=1.5, dim_trav_colum_base=1.5, tipo_concreto=1, acero_refuerzo=420000 , diametro_varilla_refuerzo = 20, varilla_columna_inferior = 3, varilla_columna_superior = 3):

        self.nombre = nombre
        self.numero_dimensiones = numero_dimensiones
        self.numero_dots = numero_dots
        self.altura_columna = altura_columna
        self.dim_travs_colum_altura = dim_travs_colum_altura
        self.dim_trav_colum_base = dim_trav_colum_base
        self.fuerza_aplicada = fuerza_aplicada
        self.resistencia_comprension = resistencia_comprension
        self.tipo_concreto = tipo_concreto
        self.acero_refuerzo = acero_refuerzo        
        self.diametro_varilla_refuerzo = diametro_varilla_refuerzo
        self.varilla_columna_inferior =  varilla_columna_inferior
        self.varilla_columna_superior = varilla_columna_superior


        ops.wipe()
        # Creación del modelo
        ops.model(self.nombre, "-ndm", self.numero_dimensiones,
                "-ndf", self.numero_dots)

    def construir_modelo(self):
        
        peso_ele = self.PES_ESP_HOR * self.dim_travs_colum_altura * \
            self.dim_trav_colum_base * self.altura_columna+self.CARGA_DISTRIBUIDA * self.altura_columna

        peso_estru = peso_ele+self.fuerza_aplicada
        
        masa = peso_estru/self.GRAVEDAD

        # Definir los materiales para la columna
        # ------------------------------------------
        # CONCRET0

        

        if self.tipo_concreto == 1:

            # Concreto inconfinado
            fpc = -self.resistencia_comprension
            epsc0 = -0.003
            fpcu = 0.2*self.resistencia_comprension
            epsu = -0.006

            ops.uniaxialMaterial("Concrete01", 1, -self.resistencia_comprension, epsc0, fpcu, epsu)

            # Concreto confinado

            # Columnas

            fcc = 1.35*fpc
            epscc0 = -0.008
            fpccu = 0.8*self.resistencia_comprension
            epsc_u = -0.020

            ops.uniaxialMaterial("Concrete01", 11, fcc, epscc0, fpccu, epsc_u)

        elif self.tipo_concreto == 2:
            # Concreto inconfinado
            fpc = -self.resistencia_comprension
            epsc0 = -0.003
            fpcu = 0.2*self.resistencia_comprension
            epsu = -0.006
            lambdac = 0.5
            ft = (0.62*(self.resistencia_comprension)**0.5)*1000
            ets = 0.1*(2*(self.resistencia_comprension)/epsc0)

            ops.uniaxialMaterial("Concrete02", 1, -self.resistencia_comprension, epsc0,
                                fpcu, epsu, lambdac, ft, ets)

            # Concreto confinado

            # Columnas

            fcc = 1.35*fpc
            epscc0 = -0.003
            fpccu = 0.2*self.resistencia_comprension
            epsc_u = -0.006

            ops.uniaxialMaterial("Concrete02", 11, fcc,
                                epscc0, fpccu, epsc_u, lambdac, ft, ets)

        elif self.tipo_concreto == 3:

            # Concreto inconfinado
            fpc = -self.resistencia_comprension
            epsc0 = -0.003
            fpcu = 0.2*self.resistencia_comprension
            epsu = -0.006
            end_strain_sitc = 0.03

            ops.uniaxialMaterial("Concrete01WithSITC",                                
                                1, -self.resistencia_comprension, epsc0, fpcu, epsu, end_strain_sitc)

            # Concreto confinado

            # Columnas

            fcc = -1.35*fpc
            epsc0 = -0.003
            fpcu = 0.2*self.resistencia_comprension
            epsu = -0.006
            end_strain_sitc = 0.03

            ops.uniaxialMaterial("Concrete01WithSITC", 11,
                                fcc, epsc0, fpcu, epsu, end_strain_sitc)

        elif self.tipo_concreto == 4:
            # Concreto inconfinado
            fpc = -self.resistencia_comprension
            epsc0 = -0.003
            fpcu = 0.2*self.resistencia_comprension
            ec_1 = (8500*(self.resistencia_comprension/1000)**0.5+110000)*1000
            beta = 0.1

            ops.uniaxialMaterial("Concrete04", 1, -self.resistencia_comprension, epsc0, fpcu, ec_1, beta)

            # Concreto confinado

            # Columnas

            fcc = 1.35*fpc
            epsc0 = -0.003
            fpcu = 0.2*self.resistencia_comprension
            ec_1 = (8500*(self.resistencia_comprension/1000)**0.5+110000)*1000
            beta = 0.1

            ops.uniaxialMaterial("Concrete04", 11, fcc, epsc0, fpcu,  ec_1, beta)

        elif self.tipo_concreto == 5:
            
            # Concreto inconfinado            
            e0 = -0.003
            n = 2
            k = 1
            alpha1 = 0.32
            ft = (0.62*(self.resistencia_comprension)**0.5)*1000
            ecr = 0.1*(2*(self.resistencia_comprension)/e0)
            b = 4
            alpha2 = 0.08

            ops.uniaxialMaterial("Concrete05", 1, -self.resistencia_comprension, e0,
                                n, k, alpha1, ft, ecr, b, alpha2)

            # Concreto confinado
            fcc = -1.35*self.resistencia_comprension
            e0 = -0.003
            n = 2
            k = 1
            alpha1 = 0.32
            ft = (0.62*(self.resistencia_comprension)**0.5)*1000
            ecr = 0.1*(2*(self.resistencia_comprension)/e0)
            b = 4
            alpha2 = 0.08

            ops.uniaxialMaterial("Concrete05", 11, -fcc,
                                e0, n, k, alpha1, ft, ecr, b, alpha2)

        elif self.tipo_concreto == 6:
            # Concreto inconfinado            
            ec = (((self.resistencia_comprension/1000)**0.25)/28)*1000
            ec_1 = (8200*(self.resistencia_comprension/1000)**(3/8))*1000
            ft = (0.62*(self.resistencia_comprension)**0.5)*1000
            et = (2*ft)/ec_1
            xp = 2
            xn = 2.3
            r = (((self.resistencia_comprension/1000)/5.2)-1.9)*1000

            ops.uniaxialMaterial("Concrete06", 1, -self.resistencia_comprension, ec,
                                ec_1, ft, et, xp, xn, r)

            # Concreto confinado
            fcc = -1.35*self.resistencia_comprension
            ec = (((self.resistencia_comprension/1000)**0.25)/28)*1000
            ec_1 = (8200*(self.resistencia_comprension/1000)**(3/8))*1000
            ft = (0.62*(self.resistencia_comprension)**0.5)*1000
            et = (2*ft)/ec_1
            xp = 2
            xn = 2.3
            r = (((self.resistencia_comprension/1000)/5.2)-1.9)*1000

            ops.uniaxialMaterial("Concrete06", 11, -fcc,
                                ec, ec_1, ft, et, xp, xn, r)
        else:
            print("No se encontro el concreto a utilizar")

        # ACERO DE REFUERZO                
        E0 = 200000000
        b = 0.007

        ops.uniaxialMaterial('Steel02', 2, self.acero_refuerzo, E0, b)

        # COLUMNAS
        # Definir la sección a la que se va a utilizar

        GJ = 1e14
        ops.section('Fiber', 2, '-GJ', GJ)

        # Definir las dimensiones

        B = self.dim_trav_colum_base
        H = self.dim_travs_colum_altura        
        cover = 0.04

        # Seccion de concreto
        self.concrete_section(11, 1, B, H, cover, 20, 5)

        # Definir el acero en la columna
        
        area_fiber = ((np.pi*self.diametro_varilla_refuerzo**2)/4)/1000**2

        # Aero inferior
                
        inicio = [-H/2+cover, -B/2+cover]
        fin = [-H/2+cover, B/2-cover]

        ops.layer('straight', 2, self.varilla_columna_inferior, area_fiber, * inicio, * fin)

        # Acero superior                
        inicio = [H/2-cover, -B/2+cover]
        fin = [H/2-cover, B/2-cover]

        ops.layer('straight', 2, self.varilla_columna_superior, area_fiber, * inicio, * fin)

        # Acero derecha
        # int(input("Ingrese el numero de varilla de la columna parte derecha: "))
        num_der = 3
        inicio = [H/2-0.20, -B/2+cover]
        fin = [H/2+0.2, B/2-cover]

        ops.layer('straight', 2, num_der, area_fiber, * inicio, * fin)

        # Acero izquierda
        # int(input("Ingrese el numero de varilla de la columna parte izquierda: "))
        num_izq = 3
        inicio = [H/2-0.20, -B/2-cover]
        fin = [H/2+0.20, B/2-cover]

        ops.layer('straight', 2, num_izq, area_fiber, * inicio, * fin)

        # Definir los nodos de la estructura

        ops.node(1, 0, 0)
        ops.node(2, 0, self.altura_columna)

        # Restricciones
        ops.fixY(0, 1, 1, 1)

        # Definir la transformacion de coordenadas
        ops.geomTransf("Linear", 1)

        # Definir los elementos de la estructura

        ops.beamIntegration(
            'HingeRadau', 2, 2, self.dim_travs_colum_altura/2, 2, self.dim_travs_colum_altura/2, 2)
        ops.element('forceBeamColumn', 1, 1, 2, 1, 2)

        # Definir los recorders
        ops.recorder("Node", "-file", self.FILE_DESPLAZAMIENTO,
                    "-time", "-node", 2, "-dof", 1, 2, 3, "disp")
        ops.recorder("Node", "-file", self.FILE_REAC_BASE, "-time",
                    "-node", 1, "-dof", 1, 2, 3, "reaction")
        
        # Definir los patrones de carga

        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)

        # Cargas

        ops.load(2, 0, (-1*masa*self.GRAVEDAD), 0)

        ops.mass(2, (masa/self.GRAVEDAD+8), 1e-8, 1e-8)

        # Crear el analisis

        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("BandGeneral")
        ops.test("EnergyIncr", 1.0e-6, 10)
        ops.algorithm("Newton")
        numero_pasos = 10
        ops.integrator("LoadControl", 1/numero_pasos)
        ops.analysis("Static")
        ok = ops.analyze(numero_pasos)

        if ok == 0:
            print("CARGA AXIAL APLICADA")
        else:
            print("ERROR AL APLICAR LA CARGA AXIAL")

        ops.loadConst("-time", 0.0)

        # Analisis modal
        num_modos = 1
        omega_sq = ops.eigen("-fullGenLapack", num_modos)
        omega_sq = np.array(omega_sq)
        omega = omega_sq**0.5
        periodo = (2*np.pi)/omega

        for i in range(num_modos):
            print("Modo", i+1, "T= ", periodo[i], "s")

        ops.remove("recorders")

        nodo = [1, 2]
        dof = [1, 2, 3]

        nodo1 = nodo[0]
        nodo2 = nodo[1]

        gra_liber1 = dof[0]
        gra_liber2 = dof[1]
        gra_liber3 = dof[2]

        # Resultados
        RESULTADOS = []

        RESULTADOS.append(str("El desplazamiento en el nodo ")+ str(nodo2) +" para el grado de libertad "+str(gra_liber1)+" es:"+str(ops.nodeDisp(2, 1))+" m")
        RESULTADOS.append(str("El desplazamiento en el nodo ")+ str(nodo2) +" para el grado de libertad "+str(gra_liber2)+" es:"+str(ops.nodeDisp(2, 2))+" m")
        RESULTADOS.append(str("El desplazamiento en el nodo ")+ str(nodo2) +" para el grado de libertad "+str(gra_liber3)+" es:"+str(ops.nodeDisp(2, 3))+" rad/s")

        RESULTADOS.append(str("La reaccción en el nodo ")+ str(nodo1) +" para el grado de libertad "+str(gra_liber1)+" es:"+str(ops.nodeReaction(1, 1))+" kN")
        RESULTADOS.append(str("La reaccción en el nodo ")+ str(nodo1) +" para el grado de libertad "+str(gra_liber2)+" es:"+str(ops.nodeReaction(1, 2))+" kN")
        RESULTADOS.append(str("La reaccción en el nodo ")+ str(nodo1) +" para el grado de libertad "+str(gra_liber3)+" es:"+str(ops.nodeReaction(1, 3))+" rad/s")
                    
        datos_1=np.loadtxt(self.FILE_DESPLAZAMIENTO)        
        datos_2=np.loadtxt(self.FILE_REAC_BASE)        
                
        serie_total = [ [ datos_1, [1,2,3],'DESPLAZAMIENTO' ], [ datos_2, [1,2,3],'REACCION_BASE' ] ]

        if serie_total:
            self.construye_figura(serie_total=serie_total, recoder='*')

        self.mostrar_resultados(lista_mensajes=RESULTADOS)


        
        
    def concrete_section(self, core_tag, cover_tag, B, H, cover, nfcore, nfcover):
        ops.patch("rect", core_tag, nfcore, 1, -H/2 +
                cover, -B/2+cover, H/2+cover, B/2+cover)
        ops.patch("rect", cover_tag, nfcover,
                1, -H/2, -B/2, -H/2+cover, B/2)
        ops.patch("rect", cover_tag, nfcover, 1, H/2-cover, -B/2, H/2, B/2)
        ops.patch("rect", cover_tag, nfcore, 1, -H /
                2+cover, -B/2, H/2-cover, B/2+cover)
        ops.patch("rect", cover_tag, nfcore, 1, -H /
                2+cover, B/2-cover, H/2-cover, B/2)


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
                plt.xlabel('Desplazamiento H. (mm)')
                plt.ylabel('Carga Aplicada (Kn)')
                print ("Normal :",i)
                plt.plot(array_mek[:,0], array_mek[:,i],label="Resultado")
                plt.legend()                           
        else :
            
            from PIL import Image
            import matplotlib.image as mpimg
            import os
            from pathlib import Path
            
            BASE_DIR = Path(__file__).resolve().parent.parent

            IMAGEN_COLUMNA = os.path.join(os.path.join(BASE_DIR,'modelos'), 'Columna_con_carga_distribuida.jpeg')

            pil_img = Image.open(IMAGEN_COLUMNA)

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
                    
                    ax3.plot(array_mek[:,0], array_mek[:,0],'tab:orange')
                    ax3.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=8)
                    ax3.set_xlabel('Desplazamiento H. (mm)')
                    ax3.set_ylabel('Carga Aplicada (Kn)')
                                        
                    ax4.plot(array_mek[:,0], array_mek[:,1],'tab:orange')
                    ax4.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}", fontsize=8)
                    ax4.set_xlabel('Desplazamiento H. (mm)')
                    ax4.set_ylabel('Carga Aplicada (Kn)')
                    
                    ax5.plot(array_mek[:,0], array_mek[:,2],'tab:orange')
                    ax5.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=8)
                    ax5.set_xlabel('Desplazamiento H. (mm)')
                    ax5.set_ylabel('Carga Aplicada (Kn)')

                elif calculo_var == 'DESPLAZAMIENTO':

                    ax6.plot(array_mek[:,0], array_mek[:,0],'tab:red')
                    ax6.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=8)
                    ax6.set_xlabel('Desplazamiento H. (mm)')
                    ax6.set_ylabel('Carga Aplicada (Kn)')
                    
                    ax7.plot(array_mek[:,0], array_mek[:,1],'tab:red')
                    ax7.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}\n", fontsize=8)
                    ax7.set_xlabel('Desplazamiento H. (mm)')
                    ax7.set_ylabel('Carga Aplicada (Kn)')

                    ax8.plot(array_mek[:,0], array_mek[:,2],'tab:red')
                    ax8.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=8)
                    ax8.set_xlabel('Desplazamiento H. (mm)')
                    ax8.set_ylabel('Carga Aplicada (Kn)')
                    
                else:            
                    pass
            
            fig.suptitle("MODELOS MATEMÁTICOS: 4.2	COLUMNA DE HORMIGÓN ARMADO CON CARGA DISTRIBUIDA A LO LARGO DE LA COLUMNA")                             
            fig.subplots_adjust(left=0.05, bottom=0.076, right=0.850, top=0.88, wspace=0.562, hspace=1)
    
        mng = plt.get_current_fig_manager()
        mng.set_window_title("MODELOS MATEMÁTICOS UTPL")
        mng.resize(*mng.window.maxsize())                
        plt.show() 

    def mostrar_resultados(self, lista_mensajes=[],ax=None):
                
        fig, ax1 = plt.subplots(1, 1, constrained_layout=True, sharey=True)       
        ax1.set_xlabel('Desplazamiento H. (mm)')
        ax1.set_ylabel('Carga Aplicada (Kn)')   
                        
        resultado = ''
        
        for index, cadena in enumerate(lista_mensajes):           
            resultado += cadena+'\n'
                        
        ax1.axes.xaxis.set_visible(False)
        ax1.axes.yaxis.set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.text(0.038, 0.75, resultado , horizontalalignment='left', verticalalignment='center_baseline', transform=ax1.transAxes, bbox={'facecolor': 'orange', 'alpha': 0.25, 'pad': 5}, fontsize=8)        
        fig.suptitle('\n\nResultados', fontsize=16)

        plt.show()    
