/* $TOG: Lower.c /main/9 1998/02/06 15:44:31 kaleb $ */

/* 
 
Copyright 1988, 1998  The Open Group

All Rights Reserved.

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
OPEN GROUP BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Except as contained in this notice, the name of The Open Group shall not be
used in advertising or otherwise to promote the sale, use or other dealings
in this Software without prior written authorization from The Open Group.

*/

/* $XFree86: xc/lib/Xmu/Lower.c,v 1.9 1999/07/11 08:49:17 dawes Exp $ */

#define  XK_LATIN1
#include <X11/keysymdef.h>
#include "CharSet.h"
#include "SysUtil.h"

#include <stdio.h>

#if 0
#ifndef HAS_SNPRINTF
#undef SCOPE
#define SCOPE static
#include "snprintf.c"
#endif
#endif

#include <stdarg.h>
#define Va_start(a,b) va_start(a,b)

/*
 * ISO Latin-1 case conversion routine
 */
#define XmuTolower(c)							 \
((c) >= XK_a && (c) <= XK_z ?						 \
 (c) : (c) >= XK_A && (c) <= XK_Z ?					 \
 (c) + (XK_a - XK_A) : (c) >= XK_Agrave && (c) <= XK_Odiaeresis ?	 \
 (c) + (XK_agrave - XK_Agrave) : (c) >= XK_Ooblique && (c) <= XK_Thorn ? \
 (c) + (XK_oslash - XK_Ooblique) :					 \
 (c))

#define XmuToupper(c)							 \
((c) >= XK_A && (c) <= XK_Z ?						 \
 (c) : (c) >= XK_a && (c) <= XK_z ?					 \
 (c) - (XK_a - XK_A) : (c) >= XK_agrave && (c) <= XK_odiaeresis ?	 \
 (c) - (XK_agrave - XK_Agrave) : (c) >= XK_oslash && (c) <= XK_thorn ?	 \
 (c) - (XK_oslash - XK_Ooblique) :					 \
 (c))

/*
 * Implementation
 */
void
XmuCopyISOLatin1Lowered(char *dst, _Xconst char *src)
{
  register unsigned char *dest, *source;

  for (dest = (unsigned char *)dst, source = (unsigned char *)src;
       *source;
       source++, dest++)
    *dest = XmuTolower(*source);
  *dest = '\0';
}

void
XmuCopyISOLatin1Uppered(char *dst, _Xconst char *src)
{
  register unsigned char *dest, *source;

  for (dest = (unsigned char *)dst, source = (unsigned char *)src;
       *source;
       source++, dest++)
    *dest = XmuToupper(*source);
  *dest = '\0';
}

int
XmuCompareISOLatin1(_Xconst char *first, _Xconst char *second)
{
  register unsigned char *ap, *bp;

  for (ap = (unsigned char *)first, bp = (unsigned char *)second;
       *ap && *bp && XmuTolower(*ap) == XmuTolower(*bp);
       ap++, bp++)
    ;

  return ((int)XmuTolower(*ap) - (int)XmuTolower(*bp));
}

void
XmuNCopyISOLatin1Lowered(char *dst, _Xconst char *src, register int size)
{
  register unsigned char *dest, *source;

  if (size > 0)
    {
      for (dest = (unsigned char *)dst, source = (unsigned char *)src;
	   *source && size > 1;
	   source++, dest++, size--)
	*dest = XmuTolower(*source);
      *dest = '\0';
    }
}

void
XmuNCopyISOLatin1Uppered(char *dst, _Xconst char *src, register int size)
{
  register unsigned char *dest, *source;

  if (size > 0)
    {
      for (dest = (unsigned char *)dst, source = (unsigned char *)src;
	   *source && size > 1;
	   source++, dest++, size--)
	*dest = XmuToupper(*source);
      *dest = '\0';
    }
}

int
XmuSnprintf(char *str, int size, _Xconst char *fmt, ...)
{
  va_list ap;
  int retval;

  if (size <= 0)
    return (size);

  Va_start(ap, fmt);

#if 0
  retval = vsprintf(str, fmt, ap);
  if (retval >= size)
    {
      fprintf(stderr, "WARNING: buffer overflow detected!\n");
      fflush(stderr);
      abort();
    }
#else
  retval = vsnprintf(str, size, fmt, ap);
#endif

  va_end(ap);

  return (retval);
}
