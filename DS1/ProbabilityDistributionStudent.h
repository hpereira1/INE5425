
#ifndef PROBABILITYDISTRIBUTIONSTUDENT_H
#define PROBABILITYDISTRIBUTIONSTUDENT_H

#include "ProbabilityDistributionBase.h"
#include "PDProxy.h"
#include "SolverProxy.h"

/*!
* Essa classe deve encontrar o valor x de uma distribuição de probabildiade cuja probabilidade acumulada corresponda a um valor dado.
* Vamos inicialmente ver o caso comum, que é achar o valor da probabildiade acumulada de uma distribuição até o valor x dado. 
* Por exemplo, para encontrar o valor da probabilidade acumulada até o valor x=1.5 numa dsitrbuição de probabilidade normal padrão (com média 0 e desvio-padrão 1),
   podemos simplesmente calcular a integral definida de -infinito até 1.5. Como "-infinito" não existe numericamente, poderíamos simplesmente colocar o menor número real representável.
* Porém, podemos tirar proveito de conhecimentos da disciplina de probabilidade e estatística e saber que praticamente 100% da área de uma distribuição de probabilidade está entre
   a média - 4xdesvio e média + 4xdesvio. Portanto poderíamos apenas calcular a integral entre -4 e 1.5. Melhor ainda seria tirar proveito do conhecimento que a probabilidade acumulada
   de uma distribuição normal até a média é 0.5 (50%) e então calcular a integral apenas de 0 a 1.5 e somar 0.5 (que é o valor pré-conhecido da integral até 0). Esse tipo de abordagem
   evita *muitos* cálculos desnecessários e pode reduzir consideravelmente o tempo de simulação de um modelo.
* Agora valor para o caso citado inicialmente: encontrar o valor x de uma distribuição de probabildiade cuja probabilidade acumulada corresponda a um valor dado.
* Por exemplo: Qual é o valor de x de uma curva normal padrão (média=0, desvio=1) tal que a probabilidade acumulada até x seja 0.4?
* Assumindo que podemos usar apenas o método "integrate" da clase "Solver" para calcular as integrais, teríamos que ficar "tentando" achar o valor "max" até que o valor da integral
   seja o valor dado. Nesse caso, max seria o valor x procurado, Isso e muito ineficiente.
* Uma abordagem mais interessante é utilizar o conhecimento disponível para definir limites iniciais aceitáveis, e então utilizar métodos numéricos análogos ao de "busca de raízes reais"
   para encontrar o valor x. No caso deste exemplo (Qual é o valor de x de uma curva normal padrão (média=0, desvio=1) tal que a probabilidade acumulada até x seja 0.4), poderíamos
   começar definindo o intervalo de busca entre -4 e 0, já que sabemos que a probabilidade acumulada nesses pois pontos é 0 e 0.5, respectivamente. Portanto o valor x cuja integral
   acumulada é 0.4 deve obrigatoriamente estar entre -4 e 0. Sendo g(x) a integral acumulada, então g(-4)=0 e g(0)=0.5. A pergunta é g(x)=0.4, qual é x?. 
* O método da da secante para busca de raízes pode ser adaptado para resolver esse problema. Ao invés de procurar x tal que f(x)=0 (se é 0, é a raiz da função), procuramos x 
   tal que g(x)=probabilidade acumulada desejada. Trata-se de um simples problema de equação de reta. Só isso. Se uma reta no ponto x=-4 passa por y=0 g(-4)=0 
   e no ponto x=0 passa pelo ponto y=0.5 g(0)=05, então em qual ponto x ela pasa por y=0.4? Basta resolver a equação da reta e ter a próxima estimativa de onde procurar pela "raiz"
   dessa função. Se o valor encontrado for o procurado (ou muito próximo dele), terminamos. Senão, definimos o novo intervalo (que é menor que o anterior) e continuamos *recursivamente* 
   a busca até encontrar a solução ou alcançar a quantidade máxima de recursões. 
*/
class ProbabilityDistributionStudent : public ProbabilityDistributionBase {
public:

    /*!dd
    * Todos os métodos "inverseXYZ" abaixo devem retornar o valor de x tal que a probabilidade acumulada da distribuição XYZ até x seja "cumulativeProbability".
    * Esses métodos devem tirar proveito do conhecimento disponível sobre cada uma das distribuições de probabilidade para definir os intervalor iniciais para a "busca de raízes", ou
    *   melhor, a *busca pela probabilidade acumulada*. 
    * Uma vez definidos esses intervalos, cada um desses métodos deve invocar uma função *recursiva* que implementa um método
    *   análogo ao da secante para busca de raízes (como explicado acima) e que refina essa busca em intervalos cada vez menores até encontrar o valor desejado ou alcançar a quantidade
    *   máxima de recursoes.
    * Esses métodos *recursivos* estão mais abaixo, e têm o nome "findInverseXYZ". Porém, os métodos "inverseXYZ" NÃO DEVEM invocá-los diretamente.
    * Para invocar um método "findInverseXYZ", deve-se usar "PDProxy::findInverseXYZ". A primeira chamada desse método recursivo deve usar o parâmetros "recursions"=0.  
    * Outro aspecto muito importante nessa impementação, é que é comum durante uma simulação que o cada um desses métodos seja invocado milhares ou milhões de vezes com eatamente os mesmos
        parâmetros. Se os parâmetros de entrada são os mesmos de uma execução anterior, então o resultado já havia sido calculado antes. Para evitar o recálculo da mesma
        coisa milhões de vezes, é imprescindível que hava uma espécie de "cache", que mapeie o conjunto de parâmetros de entrada para o valor já calculado. Assim, antes de fazer qualquer
        cálculo, esses métodos deveriam primeiro verificar nessa "cache" se o que está sendo pedido já foi previamente calculado e, se foi, simplesmente retorne o valor da cache. Se não foi,
        então o cálculo descrito deve ser realziado e o resultado deve ser guardado na cache antes de ser retornado.
    */
	static double inverseChi2(double cumulativeProbability, double degreeFreedom){
	    // DESENVOLVER
	    return 0.0;
	}
	static double inverseFFisherSnedecor(double cumulativeProbability, double d1, double d2){
	    // DESENVOLVER
	    return 0.0;
	}
	static double inverseNormal(double cumulativeProbability, double mean, double stddev){
	    // DESENVOLVER
        double invnormalresult;
        invnormalresult = cumulativeProbability*stddev + mean;
	    return invnormalresult;	    
	}
	static double inverseTStudent(double cumulativeProbability, double mean, double stddev, double degreeFreedom){
	    // DESENVOLVER
	    return 0.0;	    
	}


    /*!
    * Os métodos "findInverseXYZ" abaixo são métodos *recursivos* recebem os limites inferior "a" e superior "b" do intervalo a pesquisar,  e também o valor da 
       função (probabilidade acumulada) nesses dois pontos, que são f(a)="fa" e f(b)="fb", a quantidade de recursões "recursions" atuais, a probabilidade acumulada
       "cumulativeProbability" que se está procurando, e os demais parâmetros da distribuição de probabilidade de interesse.
    * Esses métodos devem, então, com base no método análogo ao da secante, e os valores conhecidos da probabilidade acumulada nos limites do intervalo ("fa" e "fb")
       calcular um ponto x denro desse intervalo (com base na equação de reta) que deve corresponder à probabilidade acumulada procurada.
    * Então deve calcular a probabilidade acumulada até esse valor, de maneira inteligente. Se o valor da probabilidade acumulada nesse ponto x for suficientemente
        próximo daquela procurada OU se a quantidade máxima de recursões já foi alcançada, retornar sua melhor estimativa. A quantidade máxima de iterações é "PDProxy::getMaxRecursions()".
    * Senão, com base no valor de f(x) (a probabilidade acumulada encontrada), definir um novo intervalo (que pdoe ser de a até x ou de x até b), incrementar a quantidade
        de recursões e invocar recursivamente esse método para continuar a procura.
    * Atenção: Para invocar recursivamente o método "findInverseXYZ" você DEVE invocar "PDProxy::findInverseXYZ".    
    * Como o cálculo das integrais será invocado muitas vezes recursivamente, e a própria invocação de cada método "inverseXYZ" pode ocorrer milhares ou milhões de vezes
        em cada simulação, é imprescindível fazer uma invocação inteligente das integrais, evitando qualquer recálculo desnecessário.
    * Exemplo: Suponha que esse método seja invocado com a=10, fa=0.3 e b=20. Suponha ainda que =, conforme a equação da reta, o valor estimado para x seja 15.
    * Isso significa que devemos calcular a integral acumulada até x=15. Calcular a integral de -infinito até 15 não é aceitável, uma vez que já se sabe que a integral
         acumulada até a=10 é fa=0.3. Portanto, o que deveria ser feito é calcular a integral de 10 até 15 e somar 0.3 a ela.

    * Nos métodos abaixo, então, você vai precisar calcular a integral de uma distribuição de probabilidade. Não use diretamente a sua classe "SolverStudent" para invocar
        os métodos "integrate". Ao contráio, você DEVE usar os métodos de mesmo nome, mas da classe "SolverProxy".
    * Para invocar o "integrate" você deve passar um parâmetro com o ponteiro para a função a ser integrada, que seve ser a distribuição de probabilidade de interesse.
    * Para fazer isso, você pode passar o ponteiro para um dos seguintes métodos estáticos da class abaixo.        
        class ProbabilityDistributionBase {
    	    static double chi2(double x, double degreeFreedom);
	        static double fisherSnedecor(double x, double d1, double d2);
	        static double normal(double x, double mean, double stddev);
	        static double tStudent(double x, double mean, double stddev, double degreeFreedom);
    */
    
	static double findInverseChi2(double a, double fa, double b, double fb, unsigned int recursions, double cumulativeProbability, double degreeFreedom){
	    // DESENVOLVER RECURSIVO 
	    // ...
	    // // trechos de exemplo...
	    // if (PDProxy::getMaxRecursions() ...
	    // ...
	    // // invocando o cálculo da integral de uma distribuição normal padrão de -1.5 até 1.5
	    //SolverProxy solver;
	    //double value = solver.integrate(-1.5, 1.5, &ProbabilityDistributionBase::normal, 0.0, 1.0); 
	    // ...
	    // // chamada recursiva
	    // double value = PDProxy::findInverseChi2(a, fa, b, fb, ++recursions, cumulativeProbability, degreeFreedom);
	    //
	    return 0.0;
	}
	static double findInverseFFisherSnedecor(double a, double fa, double b, double fb, unsigned int recursions, double cumulativeProbability, double d1, double d2){
	    // DESENVOLVER RECURSIVO 
	    return 0.0;	    
	}
	static double findInverseNormal(double a, double fa, double b, double fb, unsigned int recursions, double cumulativeProbability, double mean, double stddev){
	    // DESENVOLVER RECURSIVO 
	    return 0.0;	    
	}
	static double findInverseTStudent(double a, double fa, double b, double fb, unsigned int recursions, double cumulativeProbability, double mean, double stddev, double degreeFreedom){
	    // DESENVOLVER RECURSIVO 
	    return 0.0;	    
	}
};

#endif /* PROBABILITYDISTRIBUTIONSTUDENT_H */

