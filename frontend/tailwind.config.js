/** @type {import('tailwindcss').Config} */
module.exports = {
  darkMode: ["class"],
  content: [
    './pages/**/*.{js,ts,jsx,tsx,mdx}',
    './components/**/*.{js,ts,jsx,tsx,mdx}',
    './app/**/*.{js,ts,jsx,tsx,mdx}',
  ],
  theme: {
  	container: {
  		center: true,
  		padding: '2rem',
  		screens: {
  			'2xl': '1400px'
  		}
  	},
  	extend: {
  		colors: {
  			border: 'hsl(var(--border))',
  			input: 'hsl(var(--input))',
  			ring: 'hsl(var(--ring))',
  			background: 'hsl(var(--background))',
  			foreground: 'hsl(var(--foreground))',
  			primary: {
  				DEFAULT: 'hsl(var(--primary))',
  				foreground: 'hsl(var(--primary-foreground))',
  				light: 'hsl(var(--primary-light))',
  				dark: 'hsl(var(--primary-dark))'
  			},
  			secondary: {
  				DEFAULT: 'hsl(var(--secondary))',
  				foreground: 'hsl(var(--secondary-foreground))'
  			},
  			destructive: {
  				DEFAULT: 'hsl(var(--destructive))',
  				foreground: 'hsl(var(--destructive-foreground))'
  			},
  			muted: {
  				DEFAULT: 'hsl(var(--muted))',
  				foreground: 'hsl(var(--muted-foreground))'
  			},
  			accent: {
  				DEFAULT: 'hsl(var(--accent))',
  				foreground: 'hsl(var(--accent-foreground))',
  				cyan: 'hsl(var(--accent-cyan))',
  				teal: 'hsl(var(--accent-teal))',
  				coral: 'hsl(var(--accent-coral))'
  			},
  			popover: {
  				DEFAULT: 'hsl(var(--popover))',
  				foreground: 'hsl(var(--popover-foreground))'
  			},
  			card: {
  				DEFAULT: 'hsl(var(--card))',
  				foreground: 'hsl(var(--card-foreground))'
  			},
  			status: {
  				completed: 'hsl(var(--status-completed))',
  				running: 'hsl(var(--status-running))',
  				pending: 'hsl(var(--status-pending))',
  				failed: 'hsl(var(--status-failed))'
  			},
  			/* Brand Colors - New Palette */
  			twilight: {
  				'50': '#f0f0ff',
  				'100': '#e0e0ff',
  				'200': '#c2c0ff',
  				'300': '#9490ff',
  				'400': '#6560ff',
  				'500': '#3830ff',
  				'600': '#2518ff',
  				'700': '#1810d9',
  				'800': '#110B52',
  				'900': '#0a0633',
  				DEFAULT: '#110B52',
  				primary: '#110B52',
  				light: '#2518ff',
  				dark: '#0a0633'
  			},
  			seagreen: {
  				'50': '#e6ffff',
  				'100': '#b3ffff',
  				'200': '#80ffff',
  				'300': '#4dd9d9',
  				'400': '#26c7c7',
  				'500': '#00A0A0',
  				'600': '#008080',
  				'700': '#006666',
  				'800': '#004d4d',
  				'900': '#003333',
  				DEFAULT: '#00A0A0',
  				light: '#26c7c7'
  			},
  			mint: {
  				'50': '#e6fff5',
  				'100': '#b3ffe3',
  				'200': '#80ffd1',
  				'300': '#4dffbf',
  				'400': '#26ffae',
  				'500': '#00E897',
  				'600': '#00b876',
  				'700': '#008856',
  				'800': '#007F4F',
  				'900': '#003d26',
  				DEFAULT: '#00E897',
  				light: '#4dffbf',
  				turf: '#007F4F'
  			},
  			golden: {
  				'50': '#fffef0',
  				'100': '#fffcc2',
  				'200': '#fff994',
  				'300': '#fff566',
  				'400': '#f5e038',
  				'500': '#F2D513',
  				'600': '#d4b800',
  				'700': '#a69000',
  				'800': '#786800',
  				'900': '#4a4000',
  				DEFAULT: '#F2D513',
  				light: '#fff566'
  			},
  			/* Keep legacy purple/cyan/teal for compatibility */
  			purple: {
  				'50': '#f0f0ff',
  				'100': '#e0e0ff',
  				'200': '#c2c0ff',
  				'300': '#9490ff',
  				'400': '#6560ff',
  				'500': '#3830ff',
  				'600': '#2518ff',
  				'700': '#1810d9',
  				'800': '#110B52',
  				'900': '#0a0633',
  				DEFAULT: '#110B52',
  				primary: '#110B52',
  				light: '#2518ff',
  				dark: '#0a0633'
  			},
  			cyan: {
  				'50': '#e6ffff',
  				'100': '#b3ffff',
  				'200': '#80ffff',
  				'300': '#4dd9d9',
  				'400': '#26c7c7',
  				'500': '#00A0A0',
  				'600': '#008080',
  				'700': '#006666',
  				'800': '#004d4d',
  				'900': '#003333',
  				DEFAULT: '#00A0A0',
  				light: '#26c7c7'
  			},
  			teal: {
  				'50': '#e6fff5',
  				'100': '#b3ffe3',
  				'200': '#80ffd1',
  				'300': '#4dffbf',
  				'400': '#26ffae',
  				'500': '#00E897',
  				'600': '#00b876',
  				'700': '#008856',
  				'800': '#007F4F',
  				'900': '#003d26',
  				DEFAULT: '#00E897'
  			},
  			chart: {
  				'1': 'hsl(var(--chart-1))',
  				'2': 'hsl(var(--chart-2))',
  				'3': 'hsl(var(--chart-3))',
  				'4': 'hsl(var(--chart-4))',
  				'5': 'hsl(var(--chart-5))'
  			}
  		},
  		borderRadius: {
  			lg: 'var(--radius)',
  			md: 'calc(var(--radius) - 2px)',
  			sm: 'calc(var(--radius) - 4px)',
  			xl: '1rem',
  			'2xl': '1.5rem',
  			'3xl': '2rem'
  		},
  		boxShadow: {
  			'glow-sm': '0 0 15px rgba(0, 160, 160, 0.15)',
  			glow: '0 0 30px rgba(0, 160, 160, 0.2)',
  			'glow-lg': '0 0 50px rgba(0, 232, 151, 0.25)',
  			'glow-twilight': '0 0 30px rgba(17, 11, 82, 0.25)',
  			'glow-seagreen': '0 0 30px rgba(0, 160, 160, 0.25)',
  			'glow-mint': '0 0 30px rgba(0, 232, 151, 0.3)',
  			'glow-golden': '0 0 30px rgba(242, 213, 19, 0.25)',
  			'glow-cyan': '0 0 30px rgba(0, 160, 160, 0.25)',
  			'glow-teal': '0 0 30px rgba(0, 232, 151, 0.25)',
  			card: '0 4px 6px -1px rgba(17, 11, 82, 0.08), 0 2px 4px -2px rgba(17, 11, 82, 0.06)',
  			'card-hover': '0 20px 25px -5px rgba(17, 11, 82, 0.1), 0 8px 10px -6px rgba(17, 11, 82, 0.08)',
  			'inner-glow': 'inset 0 0 20px rgba(0, 160, 160, 0.1)'
  		},
  		backdropBlur: {
  			xs: '2px'
  		},
  		keyframes: {
  			'accordion-down': {
  				from: {
  					height: '0'
  				},
  				to: {
  					height: 'var(--radix-accordion-content-height)'
  				}
  			},
  			'accordion-up': {
  				from: {
  					height: 'var(--radix-accordion-content-height)'
  				},
  				to: {
  					height: '0'
  				}
  			},
  			'fade-in': {
  				from: {
  					opacity: '0'
  				},
  				to: {
  					opacity: '1'
  				}
  			},
  			'fade-out': {
  				from: {
  					opacity: '1'
  				},
  				to: {
  					opacity: '0'
  				}
  			},
  			'slide-in-right': {
  				from: {
  					transform: 'translateX(100%)'
  				},
  				to: {
  					transform: 'translateX(0)'
  				}
  			},
  			'slide-in-left': {
  				from: {
  					transform: 'translateX(-100%)'
  				},
  				to: {
  					transform: 'translateX(0)'
  				}
  			},
  			'slide-in-up': {
  				from: {
  					transform: 'translateY(100%)'
  				},
  				to: {
  					transform: 'translateY(0)'
  				}
  			},
  			'slide-in-down': {
  				from: {
  					transform: 'translateY(-100%)'
  				},
  				to: {
  					transform: 'translateY(0)'
  				}
  			},
  			'scale-in': {
  				from: {
  					transform: 'scale(0.95)',
  					opacity: '0'
  				},
  				to: {
  					transform: 'scale(1)',
  					opacity: '1'
  				}
  			},
  			'spin-slow': {
  				from: {
  					transform: 'rotate(0deg)'
  				},
  				to: {
  					transform: 'rotate(360deg)'
  				}
  			},
  			'bounce-subtle': {
  				'0%, 100%': {
  					transform: 'translateY(0)'
  				},
  				'50%': {
  					transform: 'translateY(-5px)'
  				}
  			},
  			'pulse-glow': {
  				'0%, 100%': {
  					boxShadow: '0 0 20px rgba(0, 160, 160, 0.2)'
  				},
  				'50%': {
  					boxShadow: '0 0 40px rgba(0, 232, 151, 0.4)'
  				}
  			},
  			'shimmer': {
  				'0%': {
  					backgroundPosition: '-200% 0'
  				},
  				'100%': {
  					backgroundPosition: '200% 0'
  				}
  			},
  			'float': {
  				'0%, 100%': {
  					transform: 'translateY(0)'
  				},
  				'50%': {
  					transform: 'translateY(-10px)'
  				}
  			},
  			'count-up': {
  				from: {
  					opacity: '0',
  					transform: 'translateY(10px)'
  				},
  				to: {
  					opacity: '1',
  					transform: 'translateY(0)'
  				}
  			}
  		},
  		animation: {
  			'accordion-down': 'accordion-down 0.2s ease-out',
  			'accordion-up': 'accordion-up 0.2s ease-out',
  			'fade-in': 'fade-in 0.3s ease-out',
  			'fade-out': 'fade-out 0.3s ease-out',
  			'slide-in-right': 'slide-in-right 0.3s ease-out',
  			'slide-in-left': 'slide-in-left 0.3s ease-out',
  			'slide-in-up': 'slide-in-up 0.3s ease-out',
  			'slide-in-down': 'slide-in-down 0.3s ease-out',
  			'scale-in': 'scale-in 0.2s ease-out',
  			'spin-slow': 'spin-slow 3s linear infinite',
  			'bounce-subtle': 'bounce-subtle 2s ease-in-out infinite',
  			'pulse-glow': 'pulse-glow 2s ease-in-out infinite',
  			'shimmer': 'shimmer 2s linear infinite',
  			'float': 'float 3s ease-in-out infinite',
  			'count-up': 'count-up 0.6s ease-out'
  		},
  		transitionTimingFunction: {
  			'bounce-in': 'cubic-bezier(0.68, -0.55, 0.265, 1.55)',
  			smooth: 'cubic-bezier(0.4, 0, 0.2, 1)'
  		},
  		backgroundImage: {
  			'gradient-radial': 'radial-gradient(var(--tw-gradient-stops))',
  			'gradient-conic': 'conic-gradient(from 180deg at 50% 50%, var(--tw-gradient-stops))',
  			'gradient-primary': 'linear-gradient(135deg, #110B52 0%, #00A0A0 100%)',
  			'gradient-success': 'linear-gradient(135deg, #00E897 0%, #007F4F 100%)',
  			'gradient-info': 'linear-gradient(135deg, #00A0A0 0%, #00E897 100%)',
  			'gradient-warning': 'linear-gradient(135deg, #F2D513 0%, #007F4F 100%)',
  			'gradient-danger': 'linear-gradient(135deg, #EF4444 0%, #F87171 100%)',
  			'gradient-brand': 'linear-gradient(135deg, #110B52 0%, #00A0A0 50%, #00E897 100%)',
  			'gradient-twilight-mint': 'linear-gradient(135deg, #110B52 0%, #00E897 100%)',
  			'gradient-sea-mint': 'linear-gradient(135deg, #00A0A0 0%, #00E897 100%)',
  			'gradient-mesh': 'radial-gradient(at 40% 20%, rgba(17,11,82,0.2) 0px, transparent 50%), radial-gradient(at 80% 0%, rgba(0,160,160,0.15) 0px, transparent 50%), radial-gradient(at 0% 50%, rgba(0,232,151,0.1) 0px, transparent 50%)'
  		},
  		spacing: {
  			'18': '4.5rem',
  			'88': '22rem',
  			'112': '28rem',
  			'128': '32rem'
  		},
  		fontSize: {
  			'2xs': [
  				'0.625rem',
  				{
  					lineHeight: '0.875rem'
  				}
  			]
  		},
  		zIndex: {
  			'60': '60',
  			'70': '70',
  			'80': '80',
  			'90': '90',
  			'100': '100'
  		}
  	}
  },
  plugins: [require("tailwindcss-animate")],
}
